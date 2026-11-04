"""Opt-in challenge command, separate from the historical paper estimator CLI."""
import argparse
import json
from pathlib import Path
import math
from .workflow_validation import ImageCase,challenge_plan,import_adapter,qualify_input,run_suite,write_json
from .workflow_reporting import write_report
from .workflow_design import STRATEGIES,write_design_report


def main(argv=None):
    parser=argparse.ArgumentParser(description='Controlled digital image/codebook workflow challenges; computational consistency only')
    parser.add_argument('npz',help='FZYX images/support, F by K codebook and readout/class/blank metadata')
    parser.add_argument('--output',required=True,help='new empty evidence directory')
    parser.add_argument('--adapter',default='codedcal.workflow_bardensr:make_adapter',help='module:factory; optional native dependencies are needed only for execution')
    parser.add_argument('--config',help='adapter JSON configuration')
    parser.add_argument('--policies',default='round,frame_template,barcode_model,group_search')
    parser.add_argument('--seed',type=int,default=110126)
    parser.add_argument('--strategy',choices=STRATEGIES,default='mixed_conventional',help='four mechanism-diverse probe pairs; weighted variant is experimental')
    parser.add_argument('--include-heldout',action='store_true',help='explicit research panel of five additional known digital challenges')
    parser.add_argument('--include-fractional',action='store_true',help='explicit panel of two fractional challenges with composed oracle references')
    parser.add_argument('--qualify',action='store_true',help='native/no-effect/inverse qualification only')
    parser.add_argument('--native-policy',default='identity',help='fixed complete calibration/decoder policy for all native controls')
    parser.add_argument('--run',action='store_true',help='explicitly run qualified comparisons')
    parser.add_argument('--max-seconds',type=float,default=10800)
    args=parser.parse_args(argv)
    if not math.isfinite(args.max_seconds) or args.max_seconds<=0:parser.error('--max-seconds must be finite and positive')
    policies=[p.strip() for p in args.policies.split(',')]
    if not policies or len(set(policies))!=len(policies) or any(not p or Path(p).name!=p or p in ('.','..') for p in policies):parser.error('--policies must name unique nonempty policies without path separators')
    if not args.native_policy.strip():parser.error('--native-policy must be nonempty')
    case=ImageCase.load(args.npz);out=Path(args.output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose an empty output directory; evidence is never overwritten')
    out.mkdir(parents=True,exist_ok=True);plan=challenge_plan(case,seed=args.seed,strategy=args.strategy,include_holdout=args.include_heldout,include_fractional=args.include_fractional);write_json(out/'plan.json',plan)
    design=write_design_report(out,case,plan)
    if not design['supported']:
        print('Unsupported mechanism coverage; saved design report. No adapter executions.');return 2
    if not args.run and not args.qualify:
        print('Saved plan; no workflow executions. Use --qualify or --run explicitly.');return 0
    config=json.loads(Path(args.config).read_text()) if args.config else {}
    adapter=import_adapter(args.adapter,config);gate=qualify_input(case,adapter,native_policy=args.native_policy);write_json(out/'qualification.json',gate)
    if not gate['passed']:
        write_report(out)
        print('Input/control qualification failed; comparative challenges were not run.');return 2
    if args.run:
        try:run_suite(case,adapter,policies,out/'trials',plan=plan,deadline_seconds=args.max_seconds)
        finally:summary=write_report(out)
        if not any(r['native_workflow_qualified'] for r in summary['policies'].values()):
            print('No policy passed its native workflow gate; comparative outcomes are unavailable.');return 2
    else:write_report(out)
    print('Saved native controls'+(' and specified-challenge results.' if args.run else '.'));return 0


if __name__=='__main__':raise SystemExit(main())
