# Offline reference-gate audit queued during the frozen run

Code inspection found an adjudication ordering defect: `assess` checks candidate
support before reference count, and `run_challenge` checks reference execution/
support but not count before labelling a candidate execution failure. The stated
protocol instead requires>=10 nonblank reference calls for any molecular failure
claim. Preserve raw v1 rows and enforce that unchanged requirement in an offline
adjudication after collection; do not rerun or change any fits, challenges or
thresholds. An insufficient reference must be unsupported even if the candidate
also fails. Integer references are exact native inputs, so native-qualified
integer workflows already have sufficient native counts; fractional/native-failed
panels still need the explicit correction. Test the fix before any confirmation.

Also mark spatial-first sensitivity unavailable if its matching component limit
is exceeded, while retaining independently resolved label-aware primary matches.
Do not display an unresolved spatial-first component as a computed F1 result.
Frozen execution code stays unchanged while the worker runs. No interpreted
comparative conclusion has been made from these incomplete rows.
