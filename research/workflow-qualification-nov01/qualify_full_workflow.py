"""Native controls of one fixed FULL workflow, never fresh comparative outcomes."""
from pathlib import Path
import sys,json,hashlib,signal,time
import numpy as np
from scipy.stats import chi2
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(H),str(R/'tool/src'),str(R.parent/'bardensr-upstream.nosync')]
from full_source_adapter import make_adapter
from codedcal.workflow_validation import ImageCase,qualify_input,write_json
entry=json.loads((H/'fresh-input-fullfield-512.json').read_text());p=Path(entry['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'];case=ImageCase.load(p)
out=H/'fresh-full-workflow-qualification.json';assert not out.exists()
state=json.loads((H/'development/state.json').read_text());assert state['stage']!='running','Do not compete with the frozen primary worker'
config=dict(decoder='source_cosid',source_native_scales=entry['source_conditioning_scales']);base=make_adapter(config)
class NativeView:
 def context(self,c):return base.context(c)
 def run(self,c,policy):
  assert policy=='identity'
  index=getattr(self,'index',0);self.index=index+1;folder=H/'full-workflow-control-calls';folder.mkdir(exist_ok=True)
  path=folder/f'call-{index}.json';assert not path.exists();begin=time.monotonic()
  write_json(path,dict(stage='running',call=index,native_policy='barcode_model',image_sha256=hashlib.sha256(np.ascontiguousarray(c.images).tobytes()).hexdigest()))
  try:
   result=base.run(c,'barcode_model')
   write_json(path,dict(stage='completed',call=index,native_policy='barcode_model',catalog=result.catalog.as_dict(),corrections=result.corrections,coverage=result.coverage,diagnostics=result.diagnostics,wall_seconds=time.monotonic()-begin))
   return result
  except BaseException as e:
   write_json(path,dict(stage='failed',call=index,native_policy='barcode_model',error=type(e).__name__+': '+str(e),wall_seconds=time.monotonic()-begin));raise
class Deadline(BaseException):pass
signal.signal(signal.SIGALRM,lambda *args:(_ for _ in ()).throw(Deadline()));signal.alarm(1800);start=time.monotonic()
try:
 row=qualify_input(case,NativeView());diag=row['calls'][0].get('diagnostics',{});b=diag.get('blank_calls');n=diag.get('nonblank_calls');nb=len(case.blank_classes);ng=len(np.unique(case.class_ids))-nb
 upper=float(.5*chi2.ppf(.95,2*(b+1))*ng/nb/n) if n and b is not None and nb else None
 row.update(native_policy='barcode_model',complete_workflow_control=True,config=config,input_sha256=entry['sha256'],native_blank_proxy_upper95=upper,eligible_for_confirmation=row['passed'] and upper is not None and upper<=.1,wall_seconds=time.monotonic()-start)
 write_json(out,row);print('FULL WORKFLOW',row['passed'],n,b,upper,row['eligible_for_confirmation'],flush=True)
except BaseException as e:
 write_json(out,dict(passed=False,eligible_for_confirmation=False,completed=False,native_policy='barcode_model',config=config,input_sha256=entry['sha256'],error=type(e).__name__+': '+str(e),wall_seconds=time.monotonic()-start));raise
finally:signal.alarm(0)
