# Coupled benchmark verification report

The reduced 1D benchmark passed its reaction-vector balance checks for all six implemented reactions. The zero-reaction transport limit had relative L2 and maximum errors of 1.8605e-5 and 2.0245e-5, respectively, against the analytical reference. The prescribed-pressure heterogeneous-bar limit had a relative displacement L2 error of 1.3165e-15 and normalized equilibrium residual of 3.3284e-14; the zero-pressure displacement was exactly zero.

All 12 baseline combinations (S1–S3, Kozeny-Carman or power-law permeability, and one-way or two-way coupling) passed the frozen per-step C, S, Ca, Ba, Cl, charge, and pore-volume checks after a transport positivity error was exposed and corrected. The initial implementation's charge normalisation and silent negative-inventory clipping were corrected before final reporting.

The refinement screen comprised 30 outputs across the three scenarios, two permeability laws, and five response metrics. No Stage 1 to Stage 2 difference exceeded the predeclared 1% threshold. The maximum relative Stage 1 to Stage 2 difference was 0.4258%; the maximum coarse-to-fine refinement envelope was 2.7783%. The decision rule therefore retains one-way chemistry-to-mechanics coupling for this synthetic benchmark. This is a numerical scoping decision only, not a geological inference.

The model remains unsuitable for a site-specific interpretation until porewater inputs, kinetic rates, mineral volumes, pressure history, mechanical properties and relevant boundary conditions are independently constrained. It does not resolve aqueous activity corrections, full carbonate or sulfur speciation, permeability evolution from measured textures, fracture mechanics, or fossil-scale 3D geometry.
