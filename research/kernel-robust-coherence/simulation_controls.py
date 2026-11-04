"""Frozen synthetic controls: no empirical outcomes, no seed selection."""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.special import logsumexp
from scipy.stats import beta,t as student_t
from core import collapsed,population,coherence,witness_vector
from inference import DEFAULT_PHASES,DEFAULT_SCALES

SEED=20261029
SCENARIOS=[
    ('stable',139,1500,1.,'stable',True),
    ('boundary_null',139,1500,.96,'stable',True),
    ('motion_null',139,1500,.8,'stable',True),
    ('common_complex_change',139,1500,1.,'complex',True),
    ('heterogeneous_phase',139,1500,1.,'phase',True),
    ('heterogeneous_amplitude',139,1500,1.,'amplitude',True),
    ('excluded_within_window_change',139,1500,.64,'unstable',False),
    ('excluded_motion_error_compensation',139,1500,.8,'compensation',False),
    ('large_synthetic_stable',4096,200,1.,'stable',True),
    ('large_synthetic_complex_change',4096,200,1.,'complex',True)]


def batched_inference(x,capacity,g=.96):
    """Vectorized equivalent of the tested single-experiment inference functions."""
    reps,ng,_=x.shape;logs=[]
    for theta in DEFAULT_PHASES:
        for scale in DEFAULT_SCALES:
            w=x@witness_vector(g,theta,scale)
            u=w/(capacity*(1+g*(scale+1/scale)/2))
            assert np.max(abs(u))<=1+1e-10
            for bet in [.25,.5,1.]:
                with np.errstate(divide='ignore'):logs.append(np.log1p(bet*u).sum(1))
    loge=logsumexp(logs,axis=0)-np.log(len(logs))
    radius=np.sqrt(2*np.log(8/.05)/ng);mu=x.mean(1)/capacity
    lo=mu-radius;hi=mu+radius;support=np.all(lo[:,2:]>0,axis=1)
    distance=np.maximum(np.maximum(lo[:,:2],-hi[:,:2]),0)
    lower=np.zeros(reps)
    lower[support]=np.linalg.norm(distance[support],axis=1)/np.sqrt(hi[support,2]*hi[support,3])
    w=x@witness_vector(g,0.,1.);se=w.std(1,ddof=1)/np.sqrt(ng)
    t_call=w.mean(1)-student_t.ppf(.95,ng-1)*se>0
    return dict(e=loge>=np.log(20),rectangle=lower>g,t_approx=t_call,
                positive_reference_rectangle=support,mean_log_e=loge)


def draw(rng,reps,ng,a,kind):
    amp=np.full((2,4),.85);offset=np.zeros((2,4))
    if kind=='complex':amp[:,2:]=.7;offset[:,2:]=.65
    if kind=='phase':offset[0,2:]=.55;offset[1,2:]=-.55
    if kind=='amplitude':amp[0,2:]=.5;amp[1,2:]=1.
    if kind=='unstable':amp[:]=np.array([.25,1.,1.,.25])
    phase=rng.uniform(-np.pi,np.pi,(reps,ng,2,1))
    motion=(2*rng.integers(0,2,(reps,ng,2))-1)*np.arccos(a)
    center=np.broadcast_to(phase,(reps,ng,2,4)).copy()+offset
    center[:,:,:,2:]-=motion[:,:,:,None]
    if kind=='compensation':center[:,:,:,2:]+=motion[:,:,:,None]
    # Exact binomial mean of independent bounded phasors, conditional on motion.
    balance=(2*rng.binomial(64,.5,center.shape)-64)/64
    z=(amp+1j*np.sqrt(1-amp**2)*balance)*np.exp(1j*center)
    diagonal=collapsed(z,np.ones_like(z.real)).sum(2)
    mixed=collapsed(z.sum(2),np.full((reps,ng,4),2.))-diagonal
    # Equality of two permanent independent four-way codes occurs with p=1/4.
    same=rng.random((reps,ng))<.25
    ordinary=diagonal+np.where(same,1.,-1/3)[:,:,None]*mixed
    collision=diagonal*(~same)[:,:,None]/.75
    pa=amp[:,0]*np.exp(1j*offset[:,0]);pb=amp[:,2]*np.exp(1j*offset[:,2])
    expected_g=coherence(population(np.ones((2,4)),pa,pb,a)) if kind not in ['unstable','compensation'] else None
    return dict(ordinary_H4=(ordinary,16.),source_resolved_oracle=(diagonal,2.),
                collision_aware=(collision,2/.75)),expected_g


def interval(k,n):
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,
            float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert not args.output.exists();start=time.monotonic();rng=np.random.default_rng(SEED);out=[]
    for name,ng,reps,a,kind,valid in SCENARIOS:
        totals={};logsum={};expected_g=None
        for i in range(0,reps,32):
            assert time.monotonic()-start<900,'finite compute budget exceeded'
            methods,expected_g=draw(rng,min(32,reps-i),ng,a,kind)
            for method,(x,b) in methods.items():
                r=batched_inference(x,b)
                if method not in totals:totals[method]={k:0 for k in r if k!='mean_log_e'};logsum[method]=0.
                for k in totals[method]:totals[method][k]+=int(np.sum(r[k]))
                logsum[method]+=float(r['mean_log_e'].sum())
        info={}
        for method,counts in totals.items():
            info[method]={k:dict(calls=v,rate=v/reps,binomial_95_interval=interval(v,reps)) for k,v in counts.items()}
            info[method]['average_log_e']=logsum[method]/reps
        row=dict(scenario=name,synthetic_groups=ng,repetitions=reps,true_attenuation=a,
                 population_coherence=expected_g,within_model=valid,methods=info)
        out.append(row);print(json.dumps(row),flush=True)
    payload=dict(seed=SEED,g=.96,alpha=.05,phase_grid=list(DEFAULT_PHASES),scale_grid=list(DEFAULT_SCALES),
                 results=out,elapsed_seconds=time.monotonic()-start,
                 physical_data_used=False,original_physical_comparison_passed=False,
                 note='Synthetic controls only. t diagnostics are approximate; excluded controls are assumption violations, not null coverage.')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
