"""Sparse eliminated cellular-sheaf operator with fixed boundary observations.

See MODEL_SPECIFICATION.txt. A local scalar factor has stalk ker(A_f); dividing
its equation by ||A_f|| makes its squared residual the lifted incidence energy
minimised over the factor stalk. No concentration interpretation of free q.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsmr


@dataclass
class Operator:
    D: sparse.csr_matrix
    h: np.ndarray
    x0: np.ndarray
    free: np.ndarray
    blocks: dict
    q_slice: slice
    t_slice: slice
    sizes: dict


def operator(supply, edge_senders, edge_receivers, edge_routes, flux,
             response, exposure_mean, exposure_sd, beta):
    """One FOV/programme; BMP has two source routes but one response coordinate.

    supply shape (all cells, routes); receiver indices are local 0..nr-1.
    Only positive fixed allocation edges are included. Isolated receivers remain.
    """
    supply=np.asarray(supply,float)
    response=np.asarray(response,float)
    flux=np.asarray(flux,float)
    send=np.asarray(edge_senders,int); rec=np.asarray(edge_receivers,int)
    route=np.asarray(edge_routes,int)
    if supply.ndim!=2 or not np.isfinite(supply).all() or (supply<0).any():
        raise ValueError('Supply must be a finite nonnegative matrix')
    if any(len(a)!=len(flux) for a in (send,rec,route)):
        raise ValueError('Edge arrays have different lengths')
    ns=supply.size; nq=len(flux); nr=len(response)
    if nr==0:raise ValueError('No scorable receivers')
    if nq and ((send<0).any() or (send>=len(supply)).any() or (rec<0).any()
               or (rec>=nr).any() or (route<0).any() or (route>=supply.shape[1]).any()):
        raise ValueError('Invalid endpoint or route')
    if not np.isfinite(flux).all() or (flux<=0).any():raise ValueError('Edges require positive flux')
    if nq and (supply[send,route]<=0).any():raise ValueError('Flux without supply')
    if not np.isfinite(beta) or beta<0:raise ValueError('Gain must be finite and nonnegative')
    mu=np.broadcast_to(exposure_mean,(nr,)).astype(float)
    sd=np.broadcast_to(exposure_sd,(nr,)).astype(float)
    if not np.isfinite(response).all() or not np.isfinite(mu).all() or not np.isfinite(sd).all() or (sd<=0).any():
        raise ValueError('Invalid response/reference statistics')
    qslice=slice(ns,ns+nq); tslice=slice(ns+nq,ns+nq+nr)
    rslice=slice(tslice.stop,tslice.stop+nr)
    rr=[];cc=[];vv=[];h=[];block_rows={'allocation':[],'aggregation':[],'response':[]}
    def add(coeff,rhs,block):
        row=len(h)
        norm=np.sqrt(sum(v*v for _,v in coeff))
        if norm==0:raise ValueError('Zero local equation')
        for col,val in coeff:
            rr.append(row);cc.append(col);vv.append(val/norm)
        h.append(rhs/norm);block_rows[block].append(row)
    for k in range(nq):
        sidx=int(send[k])*supply.shape[1]+int(route[k])
        a=flux[k]/supply[send[k],route[k]]
        add([(qslice.start+k,1),(sidx,-a)],0,'allocation')
    incoming=[[] for _ in range(nr)]
    for k,j in enumerate(rec):incoming[j].append(k)
    for j in range(nr):
        add([(tslice.start+j,1)]+[(qslice.start+k,-1) for k in incoming[j]],0,'aggregation')
    for j in range(nr):
        add([(rslice.start+j,1),(tslice.start+j,-beta/sd[j])],-beta*mu[j]/sd[j],'response')
    D=sparse.coo_matrix((vv,(rr,cc)),shape=(len(h),rslice.stop)).tocsr()
    exposure=np.bincount(rec,weights=flux,minlength=nr)
    x0=np.r_[supply.ravel(),flux,exposure,response]
    return Operator(D,np.array(h),x0,np.arange(qslice.start,tslice.stop),
                    {k:np.array(v,int) for k,v in block_rows.items()},qslice,tslice,
                    dict(senders=len(supply),routes=supply.shape[1],edges=nq,receivers=nr))


def solve(op):
    boundary=op.x0.copy(); boundary[op.free]=0
    A=op.D[:,op.free]; target=op.h-op.D@boundary
    result=lsmr(A,target,atol=1e-10,btol=1e-10,maxiter=10000)
    x=boundary.copy();x[op.free]=result[0]
    residual=op.D@x-op.h
    initial=op.D@op.x0-op.h
    stationarity=float(np.linalg.norm(A.T@residual)/max(1,np.linalg.norm(A.T@target)))
    if result[1] not in (0,1,2,4,5) or stationarity>1e-6:
        raise RuntimeError(f'Least squares failed: stop={result[1]}, stationarity={stationarity}')
    e0=float(initial@initial); emin=float(residual@residual)
    if emin>e0+1e-7:raise RuntimeError('Relaxation increased energy')
    q=x[op.q_slice]
    stats={'energy_fixed':e0,'energy_relaxed':emin,'solver_iterations':int(result[2]),
           'solver_stop':int(result[1]),'stationarity':stationarity,
           'negative_adjusted_fluxes':int((q < -1e-8).sum()),
           'negative_adjusted_exposures':int((x[op.t_slice] < -1e-8).sum()),
           'adjusted_flux_count':len(q), 'free_variables':len(op.free),
           'operator_rows':op.D.shape[0],'operator_columns':op.D.shape[1]}
    for block,rows in op.blocks.items():
        stats['fixed_'+block+'_energy']=float(initial[rows]@initial[rows])
        stats['relaxed_'+block+'_energy']=float(residual[rows]@residual[rows])
    return stats
