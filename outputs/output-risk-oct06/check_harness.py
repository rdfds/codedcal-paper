"""Actual-kernel coordinate/inverse replay and retained map-bank qualification."""
import numpy as np
from scipy.ndimage import gaussian_filter
import exposed as E

def main():
    E.L.guard();E.T.guard();cases=[]
    for case in E.CASES:
        maps,labels,ordinary,bases,prov=E.bank(case);assert len(maps)==len(labels) and ordinary and bases;cases.append(dict(case=E.name(case),candidates=len(maps),conventional=len(ordinary),incumbents=len(bases)))
    records=[]
    for source in ['bar','culture','tissue']:
        if source=='bar':book,negative,names=E.R.qualify.codebook();unique=list(dict.fromkeys(names));ids=np.array([unique.index(n) for n in names]);decoders=E.R.DECODERS;depth=1
        elif source=='culture':book,negative,names,_=E.A.book();ids=np.arange(book.shape[1]);decoders=E.A.DECODERS;depth=30
        else:book,_=E.T.book();negative=None;ids=np.arange(book.shape[1]);decoders=E.A.DECODERS;depth=30
        image=np.zeros((len(book),depth,128,128),np.float32);z=0 if depth==1 else 16
        for j,(y,x) in enumerate((y,x) for y in [48,58,68,78] for x in [48,58,68,78]):image[:,z,y,x]+=book[:,j]
        image=gaussian_filter(image,(0,0 if depth==1 else .6,.7,.7));image/=image.max();support=np.ones_like(image,bool);zero=np.zeros((len(book),3),int);imposed=np.zeros_like(zero);imposed[4:,1]=2;imposed[4:,2]=4
        observed=E.T.move3d(image,imposed);alpha=E.T.move3d(support,imposed)
        for decoder in decoders:
            native,q=E.decode(image,support,zero,book,source,decoder,ids,negative);correct,cq=E.decode(observed,alpha,-imposed,book,source,decoder,ids,negative);assert native==correct and len(native['words'])>=10,(source,decoder,len(native['words']));score=E.O.catalog_distance(native,correct);assert score['identity_loss']==0 and score['composition_tv']==0;records.append(dict(source=source,decoder=decoder,calls=len(native['words']),exact_inverse_catalog=True))
            print('HARNESS CONTROL',source,decoder,len(native['words']),flush=True)
    E.write(E.HERE/'harness-controls.json',dict(passed=True,map_banks_checked=13,map_banks=cases,actual_kernel_inverse_records=records,pinned_libraries_unchanged=True,decoder_harness_sha256=E.sha(E.__file__),control_sha256=E.sha(__file__)));print('PASS map banks and all actual decoder inverse controls',flush=True)
if __name__=='__main__':main()
