# Development ledger

Before any new comparative outcomes: the initial 720-call planning cap was
revised to 1,600 to retain two actual decoder families, four strong calibration
policies and four exposed inputs rather than weakening the comparator set.
The enumerated maximum is 4 inputs × 4 policies × 2 decoders × 47 logical
executions + 96 numerical/native controls = 1,600. The three-hour local worker
limit, thresholds and qualification criteria are unchanged. No data outcome
informed this budget correction.

First numerical-control run: 4 failures / 6 passes. Three failures exposed an
interface mistake: the historical matcher accepts XY only, while the new public
API promises ZYX. The new suite now has a separate component-wise 3D matcher;
the historical matcher is unchanged. One failed fractional control incorrectly
assumed that +/-0.5 shifts necessarily cause additional inverse interpolation:
in the first-readout gauge their relative corrections are integers. The control
now uses +/-0.25 and explicitly verifies the actual composed operator. No public
comparative or confirmation outcomes had been run.

## Pre-comparison engineering updates

The first native-qualification worker failed before reading inputs (`cases`
instead of the manifest's `records` key); its complete log is retained. The
corrected worker began reproducing native catalogs, starting with BARseq3 field0
correlation: 70 calls and passing native/no-effect/integer/fractional controls.
No comparative calibration/challenge outcome was read during these repairs.

Integer challenges now have first-readout shift zero and a known zero-valued
16-pixel margin. Calibration/decoder/scoring contexts preserve their original
physical pixels. This prevents controlled transport from discarding original
pixels and makes full-array integer inverse references exact. Exact identical
input/policy calls may be cached; actual executions and runtime on cache hits
are zero/lookup time, never inherited times. Geometry is shared between decoder
backends only for exact input/codebook/configuration hashes. All logical probe
budgets remain four pairs. Fractional references retain composed interpolation.
Native object adapters still use their original conditional targets and failure
qualifications, with explicit coordinate offsets for padding.

The default public TAB responses are complete but have different MD5s from the
catalogued original files. Separate original CSV requests reproduce all three
publisher MD5s. Both representations and receipts are preserved. Before any new
pixels were inspected, the fixed native qualification patch was reduced from
512 to 256 pixels per side for bounded local memory/compute; no count or method
outcome informed that choice. This is the only fresh field/preparation planned.
