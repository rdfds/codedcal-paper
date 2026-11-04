# What the references actually measure

The historical oracle uses precision-weighted binding-event means in four
finite panels. Writing their noise-free projected positions as a,b,c,d,

    That = [((a+b-c-d)/2)^2 - ((a-b)^2+(c-d)^2)/4]/2
         = [(a-c)(b-d) + (a-d)(b-c)]/4.

Independent zero-mean errors of the four means cancel in expectation. This
requires disjoint events and appropriate conditional independence; clustering
and eligibility use the same observations and can alter that expectation.
For two persistent Gaussian window positions with independent variance T,
the expectation is T. For general time paths it is a signed increment-product
functional, not necessarily a variance. Panel means integrate the time path
with event-specific weights. Unequal round durations therefore imply unequal
physical temporal supports even when panel fractions agree.

The code statistic instead uses unweighted localization sums at q=0.10 and
0.15 nm^-1, and localization-count products. For a source with panel phase
sums z_t and counts n_t, its diagonal within term is

    n_k Re(z_i conj(z_j)),

and its cross term is the average of

    n_j Re(z_i conj(z_k)), n_i Re(z_j conj(z_k)).

These are summed over the two triples in each window and averaged over the
four projection angles, with the same factors as the cubic implementation.
Assignment-assisted evaluation of these terms on the original localization
rows removes the event-mean, event-count and Fourier-versus-quadratic changes
of target. Recomputing the code cubic on those same rows further removes the
selected-versus-full population change, but introduces an oracle-dependent
selection. This is a diagnostic, not a new operational method.

Even a localization-count weighted quadratic reference is not generally the
same target as a finite-q ratio. In the stationary Gaussian special case the
diagonal ratio weights source attenuation by count products times localization
kernel factors. If all eligible sources share the same displacement law,
those factors cancel in the ratio. If source laws vary, fitting exp(-q^2 T)
to a mixture does not produce an ordinary count-weighted average variance.
Weights can vary with frequency and with which panels supply phase factors.
Round-specific motion, panel dependence and changing kernels remain real
model limitations; a quiet D statistic does not establish their absence.

The six CD20 field-derived controls characterize the original frozen oracle
using candidate clusters before source QC, their panel counts, localization
multiplicities, within-panel noise and relative geometry. Both Gaussian and
within-panel residual resampling remove unknown observed mean motion before
adding a known persistent Gaussian window displacement. They condition on
inferred source labels; merged candidates may already contain multiple real
emitters. Gaussian motion is an imposed control law, not a validated model of
these acquisitions. Thus a field-specific response is a conditional sensitivity
analysis, not an independent physical calibration certificate.

The historical SDC response has an additional problem: its fitted intercept
includes real between-window displacement of the simulated settling path.
Subtracting that intercept targets only an injected random component. Comparing
the resulting number directly with a code statistic that retains physical
settling is not a fully matched-target validation. Retain raw, corrected and
conditional calibrated views separately; the old model did not match all
observed diagnostics or selection coverage.
