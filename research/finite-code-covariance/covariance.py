"""Exact source-conditioned covariance; research/design API, not an event estimator.

Requires physical/synthetic source contributions. Never infer these from codes
or expose this as a drop-in uncertainty estimate for source-unresolved data.
"""
import itertools
import numpy as np


def check_k(k):
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, (int, np.integer)) or k < 3:
        raise ValueError('K must be an integer >= 3')
    return int(k)


def kernel(a, b, c, k):
    k = check_k(k)
    return 1. if a == b == c else (-1./(k-1) if len({a,b,c}) == 2 else 2./((k-1)*(k-2)))


class Contrast:
    """Vector of real trilinear statistics with a shared symmetric code kernel.

    features[source,feature]; factors[term,3]; weights[term,output].
    All angles, frequencies and outputs share the same source code assignment.
    """
    def __init__(self, features, factors, weights, names=None):
        if np.iscomplexobj(features) or np.iscomplexobj(weights):
            raise ValueError("real features and output coefficients required; stack Re/Im explicitly")
        self.f = np.asarray(features, float)
        raw = np.asarray(factors)
        self.w = np.asarray(weights, float)
        if (self.f.ndim != 2 or raw.ndim != 2 or raw.shape[1] != 3
                or self.w.ndim != 2 or len(raw) != len(self.w)
                or not np.issubdtype(raw.dtype, np.integer)):
            raise ValueError('invalid feature/factor/weight dimensions or indices')
        self.idx = raw.astype(int)
        if (not np.isfinite(self.f).all() or not np.isfinite(self.w).all()
                or np.any(self.idx < 0) or np.any(self.idx >= self.f.shape[1])):
            raise ValueError('nonfinite contributions or invalid feature indices')
        self.names = names or [str(i) for i in range(self.w.shape[1])]
        if len(self.names) != self.w.shape[1]:
            raise ValueError('one name per output is required')

    def evaluate(self, codes, k=4):
        """Independent aggregate implementation for one or many code assignments."""
        k = check_k(k)
        codes = np.asarray(codes)
        single = codes.ndim == 1
        if single: codes = codes[None]
        if (codes.ndim != 2 or codes.shape[1] != len(self.f)
                or not np.issubdtype(codes.dtype, np.integer)
                or np.any(codes < 0) or np.any(codes >= k)):
            raise ValueError('one valid integer code per source and assignment required')
        agg = np.zeros((len(codes), k, self.f.shape[1]))
        for j in range(len(codes)):
            np.add.at(agg[j], codes[j], self.f)
        x,y,z = (agg[..., self.idx[:,a]] for a in range(3))
        sx,sy,sz = x.sum(1),y.sum(1),z.sum(1)
        terms = (k*k*(x*y*z).sum(1) - k*((x*y).sum(1)*sz
                 +(x*z).sum(1)*sy+(y*z).sum(1)*sx) + 2*sx*sy*sz)/((k-1)*(k-2))
        result = terms @ self.w
        return result[0] if single else result

    def direct_coefficients(self):
        """Independent explicit source-index reference, for small N only."""
        n = len(self.f)
        if n > 12: raise ValueError('direct reference bounded to N <= 12')
        def a(i,j,k):
            return np.prod(self.f[[i,j,k], self.idx], axis=1) @ self.w
        mu = sum((a(i,i,i) for i in range(n)), np.zeros(self.w.shape[1]))
        pairs = []
        for i,j in itertools.combinations(range(n),2):
            pairs.append(sum((a(*v) for v in set(itertools.permutations((i,i,j)))), np.zeros_like(mu))
                         +sum((a(*v) for v in set(itertools.permutations((i,j,j)))), np.zeros_like(mu)))
        triples = [sum((a(*v) for v in itertools.permutations(ijk)),np.zeros_like(mu))
                   for ijk in itertools.combinations(range(n),3)]
        def gram(rows):
            r=np.asarray(rows).reshape(-1,len(mu));return r.T@r
        return mu,gram(pairs),gram(triples)

    def moments(self, k=4):
        """Mean and exact pair/triple covariance matrices, O(N^2 d + N L^2).

        No N^3 output allocation. Roundoff of the inclusion-exclusion triple
        Gram is retained, reported, and not silently clipped to a PSD matrix.
        """
        k=check_k(k)
        x,y,z=(self.f[:,self.idx[:,j]] for j in range(3))
        diagonal=(x*y*z)@self.w
        mu=diagonal.sum(0)
        n,d=len(self.f),self.w.shape[1]
        repeated=np.empty((n,n,d))
        for o in range(d):
            keep=np.flatnonzero(self.w[:,o])
            xx,yy,zz=x[:,keep],y[:,keep],z[:,keep];w=self.w[keep,o]
            repeated[...,o]=((xx*yy)*w)@zz.T+((xx*zz)*w)@yy.T+((yy*zz)*w)@xx.T
        both=repeated+repeated.swapaxes(0,1)
        ii=np.arange(n);both[ii,ii]=0
        b=both.reshape(-1,d);pair=b.T@b/2
        repeated[ii,ii]=0
        r=repeated.reshape(-1,d)
        all_inner=np.zeros((len(self.idx),len(self.idx)))
        xyz=(x,y,z)
        for perm in itertools.permutations(range(3)):
            all_inner+=(x.T@xyz[perm[0]])*(y.T@xyz[perm[1]])*(z.T@xyz[perm[2]])
        triple=self.w.T@all_inner@self.w-2*r.T@r-6*diagonal.T@diagonal
        triple=(triple+triple.T)/2
        covariance=pair/(k-1)+triple/((k-1)*(k-2))
        scale=max(float(np.linalg.norm(self.w.T@all_inner@self.w,ord=2)),1e-300)
        return dict(mean=mu,pair_gram=pair,triple_gram=triple,covariance=covariance,
                    triple_roundoff_scale=scale,
                    triple_min_eigenvalue=float(np.linalg.eigvalsh(triple).min()))


def four_panel(z,n,include_imag=False,include_orientations=False):
    """Implemented C_A,R_A,C_B,R_B, angle averaged, at each q.

    z[source,panel,q,angle], n[source,panel]; no synthetic source selection.
    include_imag adds the complex moment's imaginary components for covariance
    bookkeeping; the shipped real estimator uses the real outputs only.
    """
    if np.iscomplexobj(n): raise ValueError("counts must be real")
    z,n=np.asarray(z),np.asarray(n,float)
    if (z.ndim!=4 or z.shape[1]!=4 or n.shape!=z.shape[:2]
            or not np.isfinite(z).all() or not np.isfinite(n).all() or np.any(n<0)
            or z.shape[2]<1 or z.shape[3]<1):
        raise ValueError('finite source x 4-panel x q x angle phases and nonnegative counts required')
    nq,na=z.shape[2:];features=[n[:,p] for p in range(4)];lookup={}
    for q,p,a,part in itertools.product(range(nq),range(4),range(na),range(2)):
        lookup[q,p,a,part]=len(features)
        features.append(z[:,p,q,a].real if part==0 else z[:,p,q,a].imag)
    parts=('real','imag') if include_imag else ('real',)
    kinds=('C','R','O1','O2') if include_orientations else ('C','R')
    names=[f'q{q}:{win}:{kind}:{part}' for q in range(nq) for win in ('A','B')
           for kind in kinds for part in parts]
    weights={}
    def add(out,p1,p2,count,q,coefficient):
        # h is symmetric in its arguments, so phase/phase/count factor order
        # may be canonicalized without changing the aggregate statistic.
        for a in range(na):
            rr=[(0,0,1),(1,1,1)] if out.endswith('real') else [(1,0,1),(0,1,-1)]
            for u,v,sign in rr:
                key=tuple(sorted((lookup[q,p1,a,u],lookup[q,p2,a,v],count)))
                weights.setdefault(key,np.zeros(len(names)))[names.index(out)]+=coefficient*sign/na
    triples=(((0,1,2),(0,1,3)),((2,3,0),(2,3,1)))
    for q,win in itertools.product(range(nq),range(2)):
        wname=('A','B')[win]
        for i,j,k in triples[win]:
            for part in parts:
                add(f'q{q}:{wname}:R:{part}',i,j,k,q,.5)
                add(f'q{q}:{wname}:C:{part}',i,k,j,q,.25)
                add(f'q{q}:{wname}:C:{part}',j,k,i,q,.25)
                if include_orientations:
                    add(f'q{q}:{wname}:O1:{part}',i,k,j,q,.5)
                    add(f'q{q}:{wname}:O2:{part}',j,k,i,q,.5)
    return Contrast(np.column_stack(features),np.array(list(weights),int),np.array(list(weights.values())),names)


def gradients(mean):
    """Log-attenuation gradients for [C_A,R_A,C_B,R_B], positive target only."""
    ca,ra,cb,rb=np.asarray(mean,float)
    if not np.isfinite([ca,ra,cb,rb]).all() or min(ca,ra,cb,rb)<=0:
        return None
    arithmetic=np.array([1/(ca+cb),-1/(ra+rb),1/(ca+cb),-1/(ra+rb)])
    geometric=.5*np.array([1/ca,-1/ra,1/cb,-1/rb])
    return dict(arithmetic=arithmetic,geometric=geometric,
                diagnostic=np.array([1/ca,-1/ra,-1/cb,1/rb]))


def cost_variance(a,p,q,k,cost=None):
    k=check_k(k)
    return (k if cost is None else cost)*(a+p/(k-1)+q/((k-1)*(k-2)))


def phase_log_gradient(mean):
    """Gradient of log sqrt(|r_A*r_B|) with real references.

    Input [Re CA, Im CA, Re RA, Im RA, Re CB, Im CB, Re RB, Im RB].
    Reference imaginary parts are bookkeeping only; shipped phase references
    use their real parts. A zero complex numerator or nonpositive real reference
    is unsupported, rather than made usable by an arbitrary variance formula.
    """
    a=np.asarray(mean,float)
    if a.shape!=(8,) or not np.isfinite(a).all():return None
    v=np.zeros(8)
    for start in (0,4):
        re,im,ref=a[start:start+3];s=re*re+im*im
        if ref<=0 or s==0:return None
        v[start:start+3]=.5*np.array([re/s,im/s,-1/ref])
    return v
