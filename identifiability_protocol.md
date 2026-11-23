# Prefit identifiability protocol, frozen before SVD

## 1. Admissible observations

Use the published mineral pattern as a qualitative categorical observation and concretion carbonate δ13C as the second closure. The carbonate target is the published concretion mean, −13.7‰ VPDB, with ±2.1‰ treated only as the between-sample spread used to scale this screen, not as analytical precision. The organic endmember is the measured rib mean, −31.7‰ VPDB. The measured shale carbonate value, −3.0‰ VPDB, is used only as a provisional background-carbonate proxy, not as a measured seawater-DIC endmember.

Direct observations may be comparator targets only when the forward mapping and support are explicit. HI/OI, sulphate abundance, approximate sulphate δ34S, the rib sulphate-TRIS offset, and mixing-derived fractions are dependent quantities and will not be counted as independent observations. Both conflicting Conc4 total-S entries remain excluded from sulphur calibration. No value is fitted in this analysis.

## 2. Structural gate, performed before sensitivities

The selected model stage is one-way chemistry-to-mechanics. The mineral-pattern operator uses chemical mineral fields, and the δ13C operator is a post-processing mixture. Neither observes displacement, strain, stress, or stiffness. Under one-way coupling the chemical state does not receive mechanical feedback. Therefore the observation operator has an exact zero derivative with respect to effective moduli, stiffness homogenisation parameters, Biot coefficient, prescribed pore pressure, traction, and mechanical boundary conditions. These parameters are structurally unobservable by coupling topology and are excluded from the numerical SVD. They are not classified as merely weakly sensitive.

## 3. Forward observation operators

### Mineral pattern

Use the final 1D calcite and barite inventories. The source-rich region is a **synthetic bone analogue**, not registered fossil anatomy. The mask boundary is integrated over a uniform envelope `b/L ∈ [0.20, 0.30]` around the synthetic initial discontinuity at 0.25. For each boundary, zone inventory is integrated over selected cells. Detection is smoothed on log inventory fraction with logistic width 0.5. The detection threshold fraction `q` is integrated over a log-uniform range `[0.01, 0.50]` of each phase's total model inventory. The operator returns prior-averaged detection probabilities for barite in the source-rich zone, barite in the matrix analogue, and calcite in the matrix analogue. This is a nuisance-prior envelope, not an empirically measured instrument detection limit.

### Carbon isotope closure

The post-processing equation is

`δ13C_calcite = f_org(δ13C_org + ε_SRB + ε_calcite-HCO3) + (1 − f_org) δ13C_background-calcite`.

Use `ε_calcite-HCO3 = +1.0‰` as the central calcite-bicarbonate enrichment and ±0.2‰ as its screening scale, based on the inorganic precipitation experiments of Romanek et al. (1992). Use `ε_SRB = 0‰` as a central heterotrophic source-to-DIC offset with a ±2‰ screening scale. The sign and magnitude are pathway- and organism-dependent; the culture study of Londry and Des Marais (2003) found small effects for most tested heterotrophic SRB but substantial exceptions. The background endmember remains the measured shale carbonate proxy at −3.0‰ with a ±2‰ analyst-selected screening scale. These are explicit assumptions, not calibrated isotope factors or evidence that the proxy equals Toarcian seawater DIC.

## 4. Parameter coordinates and screening scales

The active chemistry vector contains log-scaled `k_SRB`, `k_SOX`, `k_SDB`, `k_calcite`, `k_barite`, aqueous diffusivity, permeability reference, and initial Ba, sulphate, and labile-OM scale factors, plus `f_org`, `ε_SRB`, `ε_calcite-HCO3`, and background δ13C. The positive parameters use one screening unit of factor 3 for rates and diffusivity, factor 10 for permeability reference, and factor 2 for concentration scales. The isotope coordinates use 0.25 for `f_org`, 2‰ for `ε_SRB`, 0.2‰ for `ε_calcite-HCO3`, and 2‰ for background δ13C. These are analyst-selected perturbation scales, not posterior priors or geological uncertainty estimates.

Run S3 for the mechanism-on sensitivity screen under both Kozeny-Carman and power-law permeability. Keep the domain, boundary conditions, duration, and starting model constants fixed. Central finite differences use a 5% fraction of each screening scale. No optimizer, likelihood maximisation, or parameter fitting is used.

## 5. Distinguishability criterion, declared before step 4

Whiten probability outputs by a unit range and δ13C by 2.1‰. Form the dimensionless sensitivity matrix from the change in these whitened outputs per one screening-scale parameter displacement. Numerical singular directions are retained as nonzero only when `σ_i / σ_max ≥ 10⁻⁸`. A direction is called distinguishable only when its singular value is at least 1, meaning a one-screening-scale perturbation along that right-singular direction moves the whitened observation vector by at least one observation-scale unit. Rank and distinguishability are reported separately. Parameter columns are not individually labelled identifiable from their marginal derivatives.

## 6. Mask-boundary envelope and numerical-artifact check

First evaluate the prior-averaged mineral operator at 20 and 40 cells. Repeat the sensitivity matrix at both resolutions using time steps satisfying the explicit diffusion stability bound (`dt = 0.2 s` at 20 cells and `dt = 0.05 s` at 40 cells). The topology result for mechanics is analytical. Chemical rank classifications must be stable to the 20-to-40-cell change under both permeability laws; if not, mark the affected singular directions numerically unresolved rather than structural.

## 7. No-fit rule

This is a local information-content screen around the existing synthetic benchmark. No parameter is estimated, no predicted value is tuned to the paper, and no result is described as fossil-specific calibration. Report the rank, right-singular identifiable combinations, numerical null space, and an independent measurement that would add sensitivity to each deficient direction.
