"""Fold-local equal-animal ridge and genuinely joint cross-lag rank-one models."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve
MODELS=('AR','VAR','AR_beh','VAR_beh','rank1','graph','persistence')
GRID=(.1,1.,10.,100.)

@dataclass
class Moments:
    n:int
    mx:np.ndarray
    my:np.ndarray
    xx:np.ndarray
    xy:np.ndarray
    yy:np.ndarray

def moments(x,y):
    return Moments(len(x),x.mean(0),y.mean(0),x.T@x/len(x),x.T@y/len(x),y.T@y/len(x))

@dataclass
class Fold:
    mx:np.ndarray
    my:np.ndarray
    sx:np.ndarray
    sy:np.ndarray
    g:np.ndarray
    c:np.ndarray
    q:np.ndarray
    n:int

def combine(items):
    """All moments are animal means: no longer animal gets more weight."""
    mx=np.mean([a.mx for a in items],0);my=np.mean([a.my for a in items],0)
    xx=np.mean([a.xx for a in items],0)-np.outer(mx,mx)
    xy=np.mean([a.xy for a in items],0)-np.outer(mx,my)
    yy=np.mean([a.yy for a in items],0)-np.outer(my,my)
    sx=np.sqrt(np.maximum(np.diag(xx),0));sy=np.sqrt(np.maximum(np.diag(yy),0))
    sx=np.where(sx>1e-12,sx,1.);sy=np.where(sy>1e-12,sy,1.)
    return Fold(mx,my,sx,sy,xx/np.outer(sx,sx),xy/np.outer(sx,sy),yy/np.outer(sy,sy),sum(a.n for a in items))

def ridge(g,c,penalty):
    return solve((g+g.T)/2+penalty*np.eye(len(g)),c,assume_a='pos',check_finite=False)

def cross_rank1(f,h,alpha,maxiter=200,tol=1e-8):
    p=len(f.mx);rho=alpha/f.n;g=f.g[:2*h,:2*h];c=f.c[:2*h]
    full=ridge(g,c,rho);cross=np.stack([full[h:,0],full[:h,1]])
    u,s,vt=np.linalg.svd(cross,full_matrices=False);k=vt[0].copy();a=u[:,0]*s[0]
    b=np.zeros((2*h,2));prev=np.inf;converged=False;objectives=[]
    for iteration in range(1,maxiter+1):
        # Given a temporal kernel, each receiver has its own history plus one cross port.
        for j in range(2):
            own=np.arange(j*h,(j+1)*h);sender=np.arange((1-j)*h,(2-j)*h)
            transform=np.zeros((2*h,h+1));transform[own,np.arange(h)]=1;transform[sender,h]=k
            coef=ridge(transform.T@g@transform,transform.T@c[:,j],rho)
            b[own,j]=coef[:h];a[j]=coef[-1];b[sender,j]=a[j]*k
        # Given own histories/gains, minimise exact ridge objective over the shared kernel.
        kg=np.zeros((h,h));kc=np.zeros(h)
        for j in range(2):
            own=np.arange(j*h,(j+1)*h);sender=np.arange((1-j)*h,(2-j)*h)
            kg+=a[j]**2*g[np.ix_(sender,sender)]
            kc+=a[j]*(c[sender,j]-g[np.ix_(sender,own)]@b[own,j])
        if np.sum(a*a)<1e-24:
            k=np.zeros(h);converged=True
        else:
            k=ridge(kg,kc,rho*np.sum(a*a))
        norm=np.linalg.norm(k)
        if norm>1e-12:k/=norm;a*=norm
        else:k=np.zeros(h);a[:]=0
        for j in range(2):b[(1-j)*h:(2-j)*h,j]=a[j]*k
        objective=float(np.trace(b.T@g@b)-2*np.sum(b*c)+np.trace(f.q)+rho*np.sum(b*b))
        objectives.append(objective)
        if np.isfinite(prev) and abs(prev-objective)<=tol*max(abs(prev),1e-12):converged=True;break
        prev=objective
    result=np.zeros((p,2));result[:2*h]=b
    return result,{'converged':bool(converged),'iterations':iteration,'objective':objective,'maximum_objective_increase':float(max(np.diff(objectives),default=0.))}

def fit(f,h,alpha,model,graph=(True,True)):
    p=len(f.mx);b=np.zeros((p,2));info={'converged':True,'iterations':1}
    if model=='rank1':b,info=cross_rank1(f,h,alpha)
    elif model=='persistence':
        # Exact identity forecast in original units, independent of training estimates.
        b[0,0]=f.sx[0]/f.sy[0];b[h,1]=f.sx[h]/f.sy[1]
        offset=np.array([(f.mx[0]-f.my[0])/f.sy[0],(f.mx[h]-f.my[1])/f.sy[1]])
        return b,offset,info
    else:
        for j in range(2):
            own=list(range(j*h,(j+1)*h));both=list(range(2*h));beh=list(range(2*h,p))
            idx=own if model=='AR' else both if model=='VAR' else own+beh if model=='AR_beh' else both+beh if model=='VAR_beh' else both if graph[j] else own
            b[idx,j]=ridge(f.g[np.ix_(idx,idx)],f.c[idx,j],alpha/f.n)
    return b,np.zeros(2),info

def mse_from_moments(m,f,b,offset):
    """Original-unit predictions from sufficient statistics, not cancellation of raw baselines."""
    coef=(b/f.sx[:,None])*f.sy[None,:]
    intercept=f.my+f.sy*offset-f.mx@coef
    losses=np.diag(m.yy)-2*np.sum(coef*m.xy,axis=0)+np.diag(coef.T@m.xx@coef)+2*intercept*(m.mx@coef-m.my)+intercept**2
    if np.any(losses < -1e-7*np.maximum(np.diag(m.yy),1)):
        raise ValueError('Negative MSE beyond floating-point tolerance')
    losses=np.maximum(losses,0)
    return losses,losses/f.sy**2

def predict(x,f,b,offset):
    return (((x-f.mx)/f.sx)@b+offset)*f.sy+f.my
