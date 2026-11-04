"""Describe sparse blank controls; illustrative iid limit, never a true-FDR CI."""
from scipy.stats import beta
from condition import HERE,read,write,sha

classes=read(HERE/'codebook-classes.json')
nblank=sum(classes['blank']);ngene=len(classes['blank'])-nblank
assert (nblank,ngene)==(8,1762)

def limit(blanks,genes):
    n=blanks+genes
    if not n or not genes:return None
    # Exact one-sided binomial limit conditional on all calls, assuming iid
    # blank/nonblank events. Peaks can violate iid; code-specific selection can
    # violate exchangeability, so this is not molecular-FDR certification.
    q=float(beta.ppf(.95,blanks+1,genes))
    return q/(1-q)*(ngene/nblank)

records=[]
for i in [0,1]:
    for d,r in read(HERE/f'field-{i}/qualification.json')['decoders'].items():
        records.append(dict(field=i,decoder=d,nonblank=r['nonblank_calls'],blanks=r['blank_calls'],
            point_proxy=r['blank_proxy'],illustrative_iid_one_sided_95_upper_adjusted_rate_ratio=limit(r['blank_calls'],r['nonblank_calls'])))
write(HERE/'blank-uncertainty.json',dict(native_gates=records,blank_classes=nblank,nonblank_code_classes=ngene,
    scope='Illustrative iid binomial sampling limit for adjusted blank/nonblank rate ratios, not a molecular-FDR interval. Codebook exchangeability and independent peaks are not established. Ratios can exceed one.',
    unchanged_gate='Fixed operational native qualification uses count/support/point proxy; does not certify proxy <=.10 with confidence.',
    script_sha256=sha(__file__)))
for r in records:print(r,flush=True)
