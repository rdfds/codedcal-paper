"""Known corruptions of real native catalogs; never used to select a winning method."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(R/'tool/src'))
from codedcal.workflow_validation import Catalog,compare,write_json
out=H/'evaluation-controls.json';assert not out.exists();records=[]
for p in sorted((H/'native-qualification-v3').glob('*.json')):
 row=json.loads(p.read_text())
 if not row['passed']:continue
 native=Catalog.from_dict(row['calls'][0]['catalog']);n=len(native.classes);rng=np.random.default_rng(110126)
 altered=Catalog(native.points,native.classes[rng.permutation(n)]);shifted=Catalog(native.points+np.array([0,8,0]),native.classes)
 composed=Catalog(shifted.points-np.array([0,8,0]),shifted.classes)
 composition=Catalog(native.points,np.repeat(native.classes[0],n))
 a=compare(native,altered);g=compare(native,shifted);c=compare(native,composition);oracle=compare(native,composed)
 passed=a['identity_loss']>.05 and a['composition_tv']==0 and g['identity_loss']>.05 and g['composition_tv']==0 and c['composition_tv']>.05 and oracle['identity_loss']==0 and oracle['composition_tv']==0
 records.append(dict(input=p.name,native_calls=n,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),passed=passed,count_preserving_label_permutation=a,known_geometric_displacement=g,composition_corruption=c,oracle_correction=oracle))
write_json(out,dict(passed=bool(records) and all(r['passed'] for r in records),records=records,uses_actual_native_catalogs=True,additional_workflow_executions=0,selection_of_endpoint_preceded_comparative_outcomes=True))
print('EVALUATION CONTROLS',len(records),all(r['passed'] for r in records))
