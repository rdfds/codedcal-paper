"""Complete fixed-plane light gates after documented expensive-run interruption."""
import json
import traceback
import numpy as np
import adapter as A

def main():
    A.guard();codes,negative,*_=A.book();records=[]
    for k in range(5):
        source=A.CACHE/'conditioned-Position351'/f'crop-{k}.npz'
        with np.load(source) as data:images=data['images'];support=data['support']
        for decoder in A.DECODERS:
            path=A.HERE/'native'/f'crop-{k}__{decoder}.json'
            if path.exists():r=json.loads(path.read_text())
            elif decoder=='regression':r=dict(status='interrupted' if k==1 else 'not_run',crop=k,decoder=decoder,reason='sampling diagnosis; no numerical outcome inferred')
            else:
                try:
                    xy,calls,runtime=A.decode(images,codes,decoder)
                    r=dict(status='completed',crop=k,decoder=decoder,xy=xy.tolist(),words=calls.tolist(),coding=int((~negative[calls]).sum()),negative=int(negative[calls].sum()),support=float(np.all(support[:,13:16,24:104,24:104],axis=(0,1)).mean()),runtime=runtime)
                except Exception:r=dict(status='failed',crop=k,decoder=decoder,traceback=traceback.format_exc())
                r.update(input_sha256=A.sha(source),freeze_sha256=A.sha(A.HERE/'freeze.json'));A.write(path,r)
            records.append(r);print(k,decoder,r['status'],r.get('coding'),r.get('negative'),flush=True)
    groups={}
    for decoder in A.DECODERS:
        rows=[r for r in records if r['decoder']==decoder];complete=all(r['status']=='completed' for r in rows)
        n=sum(r.get('coding',0) for r in rows);b=sum(r.get('negative',0) for r in rows);proxy=b/n*998/26 if n else None
        groups[decoder]=dict(completed=complete,coding=n,negative=b,scaled_proxy=proxy,qualified=complete and n>=50 and proxy<=.1 and min(r.get('support',0) for r in rows)>=.99)
    A.write(A.HERE/'state.json',dict(stage='failed_fixed_plane_light_gates_with_partial_regression',planned=15,completed=sum(r['status']=='completed' for r in records),interrupted=sum(r['status']=='interrupted' for r in records),not_run=sum(r['status']=='not_run' for r in records),failed=sum(r['status']=='failed' for r in records),groups=groups,records=records,adapter_qualified=False,independent_comparison_completed=False))

if __name__=='__main__':main()
