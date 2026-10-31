"""Response-blind spatial allocations; this kernel method is not COMMOT."""
import importlib.util
import numpy as np
from scipy.spatial.distance import cdist


ROUTES={
    'IFNG':[('IFNG',('IFNGR1','IFNGR2'))],
    'TNF':[('TNF',('TNFRSF1A',))],
    'BMP':[('BMP2',('BMPR1A','BMPR2')),('BMP6',('ACVR1','BMPR2'))],
}
INPUT_GENES=sorted({g for routes in ROUTES.values() for ligand,receptor in routes for g in (ligand,*receptor)})


def capacities(h,genes,receiver_ids,programme):
    ix={g:i for i,g in enumerate(genes)}
    routes=ROUTES[programme]
    supply=np.column_stack([h[:,ix[l]] for l,_ in routes])
    receptor=np.column_stack([np.min(h[receiver_ids][:,[ix[g] for g in r]],axis=1) for _,r in routes])
    shared=h[receiver_ids,ix['BMPR2']] if programme=='BMP' else None
    return supply,receptor,shared


def cap_flux(q,supply,receptor,shared=None):
    """Downward-only route/column and shared-budget caps, preserving sparsity."""
    q=np.asarray(q,float).copy()
    if q.shape!=(len(supply),len(receptor),supply.shape[1]):raise ValueError('Flux dimensions differ')
    if not np.isfinite(q).all() or (q<0).any():raise ValueError('Invalid flux')
    # Numerical safety also guards optional external allocation backend.
    totals=q.sum(axis=1)
    q*=np.minimum(1,np.divide(supply,totals,out=np.ones_like(totals),where=totals>0))[:,None,:]
    totals=q.sum(axis=0)
    q*=np.minimum(1,np.divide(receptor,totals,out=np.ones_like(totals),where=totals>0))[None,:,:]
    if shared is not None:
        incoming=q.sum(axis=(0,2))
        q*=np.minimum(1,np.divide(shared,incoming,out=np.ones_like(incoming),where=incoming>0))[None,:,None]
    return q


def allocate(h,genes,coords,receiver_ids,programme,radius,backend='kernel'):
    receiver_ids=np.asarray(receiver_ids,int)
    supply,receptor,shared=capacities(h,genes,receiver_ids,programme)
    distance=cdist(coords,coords[receiver_ids])
    # Coincident distinct centres are allowed; self identity, not distance alone, excludes autocrine.
    support=distance<=radius
    support[receiver_ids,np.arange(len(receiver_ids))]=False
    if backend=='kernel':
        weight=np.exp(-.5*(distance/(radius/2))**2)*support
        raw=weight[:,:,None]*receptor[None,:,:]
        denom=raw.sum(axis=1)
        q=raw*np.divide(supply,denom,out=np.zeros_like(supply),where=denom>0)[:,None,:]
    elif backend=='commot':
        if not all(importlib.util.find_spec(m) for m in ('commot','anndata')):
            raise RuntimeError('COMMOT backend requires commot and anndata; use explicit --allocator kernel if unavailable')
        import anndata
        import commot
        import pandas as pd
        # COMMOT computes a cell x cell distance matrix. Mask non-targets/self before OT.
        full=cdist(coords,coords)
        allowed=np.zeros_like(full,dtype=bool); allowed[:,receiver_ids]=support
        full[~allowed]=np.inf
        adata=anndata.AnnData(X=np.asarray(h,float),var=pd.DataFrame(index=genes))
        adata.obsm['spatial']=np.asarray(coords,float)
        adata.obsp['spatial_distance']=full
        routes=ROUTES[programme]
        db=pd.DataFrame([[l,'_'.join(r),programme] for l,r in routes],columns=['ligand','receptor','pathway'])
        commot.tl.spatial_communication(adata,database_name='pilot',df_ligrec=db,
            heteromeric=True,heteromeric_rule='min',dis_thr=radius,
            cot_eps_p=.1,cot_rho=10,cot_nitermax=10000,smooth=False,pathway_sum=False)
        q=np.zeros((len(h),len(receiver_ids),len(routes)))
        for k,(l,r) in enumerate(routes):
            mat=adata.obsp['commot-pilot-'+l+'-'+'_'.join(r)]
            q[:,:,k]=mat[:,receiver_ids].toarray() if hasattr(mat,'toarray') else np.asarray(mat)[:,receiver_ids]
        if (q[~support]>1e-8).any():raise RuntimeError('COMMOT violated explicit spatial/receiver mask')
        q[~support]=0
    else:raise ValueError('Unknown allocation backend')
    q=cap_flux(q,supply,receptor,shared)
    tolerance=1e-8
    assert np.max(q.sum(axis=1)-supply,initial=0)<=tolerance
    assert np.max(q.sum(axis=0)-receptor,initial=0)<=tolerance
    assert not (q[~support]>0).any()
    if shared is not None:assert np.max(q.sum(axis=(0,2))-shared,initial=0)<=tolerance
    send,rec,route=np.nonzero(q>0)
    return dict(supply=supply,send=send,rec=rec,route=route,flux=q[send,rec,route],
                exposure=q.sum(axis=(0,2)),distance=distance[send,rec],
                supported_pairs=int(support.sum()),receptor=receptor)
