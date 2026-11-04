# Fixed model-based controls, before execution

29 October 2026. Seed 20261029. No acquisition outcomes are inputs.
These are finite implementation/power controls, not physical validation.

Use independent groups, each with two sources and 64 independent phase samples
per panel. Unit-normalized counts are one in every panel; this common scale
changes no ratios. Source phases are uniform, source displacement is independently
+/-acos(a), and panel errors are independent +/-acos(kernel magnitude).
Binomial sampling generates the exact finite phase means, avoiding a normal
approximation. Panel kernels are stable within each window in valid controls.
Every source reference power remains at least .25.

Primary g=.96, alpha=.05, 139 groups, 1500 repetitions per scenario. Eight
139-group scenarios: stable a=1; boundary null a=.96; motion null a=.8;
common complex kernel change a=1; heterogeneous phase kernels a=1;
heterogeneous amplitudes a=1; excluded within-window change a=.64; and
excluded motion/error compensation a=.8. The last two deliberately violate
the theorem and must not be counted as type-I-error controls.

Two positive power controls repeat stable a=1 and common complex kernel change
with 4096 **synthetic** groups and 200 repetitions. They do not increase the
number of independent groups in the real project. Total ten scenarios.

Compare ordinary H4, source-resolved oracle and independently design-flagged
collision-aware estimator; keep all groups including zero collision-aware
contributions. Report finite fixed-mixture e-test, simultaneous Hoeffding
rectangle, and explicitly approximate fixed-direction t diagnostic separately.
The fixed direction is theta=0 and scale=1, including when a complex kernel
change makes it unhelpful. Never fit the scored direction or select code seeds.
For the mixture, repair the initial coarse four-phase grid to sixteen phases
and five log-spaced scales [.5,sqrt(.5),1,sqrt(2),2], with three fixed bets.
This resolution is frozen here; it is not selected from simulation winners.
The primary mixture is the sole nominal family for power reporting; the
rectangle and t diagnostics are additional separately labeled summaries.

Use <=32-repetition working batches, <=15 minutes, <=32MiB retained output.
Report binomial Monte Carlo uncertainty and negative/abstaining results. A low
false-call rate in simulations does not establish the premises in real images.
Report population coherence and the true a so conservatism and assumption
violations are visible. A failure of a conservative test at 139 groups cannot
be repaired by calling the t approximation an exact certificate.
