# Second-stage screening specification

The original frozen observation and sensitivity protocol remains unchanged. These additions extend it without fitting or altering the benchmark's kinetic baseline.

## Barite model-form comparisons

The acidity switch multiplies calcite precipitation kinetics by a smooth inhibitor, `1 / (1 + 10^((6.5 − pH)/0.25))`. The switch is a kinetic hypothesis, not a measured pH threshold. A local Ba source is introduced as a charge-balanced BaCl2 inventory input to the synthetic bone-proxy zone, with the source tracked in the Ba and Cl budgets. The original benchmark condition, uniform-Ba control, acidity-only condition, local-source condition at benchmark background, and combined condition are compared. A separate low-background-Ba stress case uses a background/inlet concentration scale of 10−3 relative to the benchmark and a source-rate sweep from 10−5 to 6×10−4 synthetic concentration units s−1. That stress case lies outside the continuous Ba screening envelope of 0.5–2 and is not treated as a likely or admissible fossil condition.

For the categorical pattern screen, bone localisation is counted only when prior-averaged detection probability is at least 0.80 in the bone-proxy mask and no more than 0.20 in the matrix-proxy mask. The detection operator retains the original mask-boundary and log-uniform threshold envelopes. Because the threshold is a fraction of total model mineral inventory, this criterion tests relative spatial pattern, not absolute mineral detectability or abundance.

## LHS rank and isotope propagation

A paired 32-point Latin-hypercube design with seed 20261111 spans each of the 15 normalized parameter coordinates uniformly from −1 to +1 screening units. This is a coverage design, not a Bayesian prior or geological probability distribution. The LHS uses the 20-cell synthetic domain, S3, and both permeability laws with the same parameter points. Central-difference step remains 0.05 screening units. For points whose largest central-difference diffusivity exceeds 2×10−9 model units, dt is reduced to 0.1 s; otherwise dt=0.2 s. One timestep is held fixed across all parameter columns within each local sensitivity matrix. The original relative numerical-rank cutoff of 10−8 and practical singular-value cutoff of 1 remain in force.

The isotope envelope uses the already recorded coordinate ranges: organic fraction 0.25–0.75, SRB substrate-to-DIC offset −2 to +2‰, calcite-bicarbonate enrichment 0.8–1.2‰, background carbonate/DIC proxy −5 to −1‰, and organic endmember −32.2 to −31.2‰. These intervals propagate the stated screening scales; the resulting range is not a posterior interval. The 2.1‰ δ13C scale remains the published between-sample spread, not analytical precision.
