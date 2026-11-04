"""Lossless typed projections of already exposed inputs; no comparative fitting."""
from pathlib import Path
import csv,hashlib,itertools,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
DEST=ROOT.parent/'workflow-qualification-nov01.nosync';DEST.mkdir(exist_ok=True)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def names_map(names,negative):
    classes=list(dict.fromkeys(names));mapping=np.array([classes.index(n) for n in names])
    return mapping,np.unique(mapping[np.array(negative,bool)]),classes

def save(name,prep,images,support,book,names,negative,metadata,sources):
    mapping,blanks,classes=names_map(names,negative);f=len(images);rounds=f//4
    p=DEST/(name+'.npz');assert not p.exists(),'Retain projected evidence'
    np.savez_compressed(p,images=np.asarray(images,np.float32),support=np.asarray(support,np.float32),codebook=np.asarray(book,np.float32),
        rounds=np.repeat(np.arange(rounds),4),colors=np.tile(np.arange(4),rounds),class_ids=mapping,blank_classes=blanks,
        name=np.array(name),preparation=np.array(prep),metadata_json=np.array(json.dumps(metadata)))
    return dict(name=name,preparation=prep,path=str(p),sha256=sha(p),shape=list(images.shape),class_names=classes,
        source_receipts=sources,metadata=metadata)

def main():
    manifest=[];bar=ROOT/'outputs/barseq3-native-repair-oct05';cached=ROOT.parent/'barseq3-inputs-oct05.nosync';classes=read(bar/'codebook-classes.json')
    book=np.zeros((40,len(classes['names'])),np.float32)
    for j,s in enumerate(classes['sequences']):
        for r,c in enumerate(s):book[4*r+'GTCA'.index(c),j]=1
    for k in [0,1]:
        images=[];supports=[];source=[]
        for r in range(10):
            p=cached/f'native-repair-source-{k}-round-{r}.npz';expected=read(bar/f'field-{k}'/f'round-{r}-demons.json')['output_sha256'];assert sha(p)==expected
            with np.load(p,allow_pickle=False) as d:images.extend(d['images']);supports.extend(d['support'])
            source.append(dict(path=str(p),sha256=expected))
        meta=dict(target='native_pixel_peaks',context=[[0,1],[16,112],[16,112]],scoring_region=[[0,1],[24,104],[24,104]],
            coverage_region=[[0,1],[24,104],[24,104]],historical_reference={dec:str(bar/'recovery'/f'field-{k}'/f'reference-{dec}.json') for dec in ['correlation','regression']})
        manifest.append(save('barseq3-'+str(k),'barseq3-adult-cb',np.array(images)[:,None],np.array(supports)[:,None],book,classes['names'],classes['blank'],meta,source))
    books=read(ROOT/'outputs/barseq-public-substitute-oct06/starfinder-color-codebooks.json')['source_books']
    for kind,index in [('culture',1),('tissue',0)]:
        src=books[index];words=[tuple(int(c)-1 for c in row['color_word']) for row in src['items']];rounds=len(words[0]);absent=[]
        if kind=='culture':absent=[w for w in itertools.product(range(4),repeat=rounds) if np.bitwise_xor.reduce(w)==0 and w not in words]
        names=[row['name'] for row in src['items']]+['no_gene_'+''.join(str(c+1) for c in w) for w in absent];negative=[False]*len(words)+[True]*len(absent)
        book=np.zeros((rounds*4,len(names)),np.float32)
        for j,w in enumerate(words+absent):book[np.arange(rounds)*4+list(w),j]=1
        if kind=='culture':
            p=ROOT.parent/'starfinder-inputs-oct06.nosync/local-trace-Position351/crop-4.npz'
            expected=read(ROOT/'outputs/starfinder-public-confirmation-oct06/quality-catalog/native/crop-4__correlation.json')['input_sha256']
        else:
            p=ROOT.parent/'starfinder-tissue-inputs-oct06.nosync/fullfield-conditioned-crop-4.npz'
            expected=read(ROOT/'outputs/starfinder-public-confirmation-oct06/fullfield-tissue/crop-4.json')['array_sha256']
        assert sha(p)==expected
        with np.load(p,allow_pickle=False) as d:images=d['images'];support=d['support']
        meta=dict(target='supported_source_object_catalog',object_adapter=kind,context=[[3,27],[16,112],[16,112]],
            scoring_region=[[4,26],[24,104],[24,104]],coverage_region=[[4,26],[24,104],[24,104]],
            conditional_object_target=True,original_whole_interface_gate_failures_preserved=True,
            model_projection=kind=='culture')
        manifest.append(save('starfinder-'+kind,'starfinder-'+kind,images,support,book,names,negative,meta,[dict(path=str(p),sha256=expected)]))
    (HERE/'exposed-inputs.json').write_text(json.dumps(dict(records=manifest,all_preparations_previously_exposed=True,independent_confirmation=False),indent=2)+'\n')
    print([(r['name'],r['shape']) for r in manifest])

if __name__=='__main__':main()
