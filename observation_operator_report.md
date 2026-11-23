# Observation operator and identifiability flags

## Definition

The observation operator maps a simulated state to the quantity and support represented by a measurement. The mapping must preserve the distinction between a measured specimen value, a value calculated from measurements, and an interpretation proposed by the source paper. A data row is not a model input merely because it is published.

The operator catalogue contains 21 observable groups. `published_observation_ledger.csv` retains the 504 sample-level rows from the source-trace register with source file, sheet, cell, supplied units and a record-level evidence flag. `source_conflict_ledger.csv` retains both conflicting Conc4 sulphur values and the other source discrepancies without overwriting either source.

## Operators available in this benchmark

The mineral operator integrates calcite or barite inventory over a supplied set of cell indices (M):

`Y_phase(M) = sum[n_phase(i) * dx(i), i in M]`,

where `n_phase` is in mol m-3 of bulk volume and `dx` is in metres, giving mol m-2. A mineral-presence probability then averages a smoothed detection threshold across a stated analyst-selected prior and the source-rich/matrix mask boundary. The boundary interval and threshold prior are explicitly recorded in `identifiability_protocol.md`. They are not measured instrument detection limits, and the resulting probabilities are not geological predictions because the published specimen is not registered to this synthetic 1D grid.

The second operator is a post-processing δ13C mixture, `δ13C_calcite = f_org(δ13C_org + ε_SRB + ε_calcite-HCO3) + (1−f_org)δ13C_background-calcite`. It uses a fixed explicit central calcite-bicarbonate enrichment and a variable SRB offset. This is an isotope closure scenario, not an isotope-resolved reactive transport calculation. The shale carbonate value is used only as an assumed background proxy, not as a direct seawater-DIC measurement.

The model's current 1 mm domain and synthetic two-part initial concentration pattern do not represent measured fossil geometry or registered bone, matrix and shale compartments. Their spatial mineral profiles can be inspected as solver outputs, but cannot be compared point-by-point with Jian et al.'s specimen.

## Direct measurements

The source ledger identifies sample measurements including Rock-Eval S1, S2, TOC and Tmax; carbonate and organic-matter carbon isotopes; carbonate oxygen isotopes; TIC; total sulphur; TRIS abundance; total sulphur isotope; TRIS sulphur isotope; total nitrogen; and total nitrogen isotope. Petrographic mineral occurrences such as barite in marrow cavities, calcite in matrix, and fluorapatite replacement are treated as qualitative direct observations. The mineral operator is implemented on synthetic model masks, and δ13C is available only through a post-processing closure. The model still lacks isotope pools, full sulphur mineral closure, nitrogen, thermal maturation, specimen registration, and quantitative phase mapping.

## Derived values and dependent observations

HI and OI are calculated Rock-Eval indices. Sulphate abundance is calculated as total sulphur minus TRIS. The reported sulphate δ34S values are approximate mass-balance calculations, not direct isotope measurements. The rib sulphate-TRIS offset depends on these derived sulphate values and paired sample treatment. Carbon-derived DIC fractions and inferred water δ18O values are also derived quantities whose values depend on selected endmembers, temperature, and fractionation assumptions. They are flagged as dependent or secondary checks and must not be entered as independent likelihood terms alongside their source measurements.

For the rib sulphur isotope comparison, the final article's TRIS standard deviation is ±2.6‰; the source workbook's ±2.0‰ remains recorded as a discrepancy. Conc4-A total sulphur is 6.3559 wt.% in the article table and 9.5 wt.% in the Rock-Eval workbook. Conc4-B is 5.6403 wt.% and 7.4535 wt.%, respectively. Both pairs remain unresolved and are excluded from total-sulphur mass-balance calibration.

## Non-unique mappings

Observed barite and calcite distributions do not uniquely identify saturation, local supply, nucleation kinetics, growth rates, or later alteration. Carbonate δ13C does not uniquely identify organic-carbon fraction without source endmembers and microbial/mineral fractionations. Carbonate δ18O does not independently separate water composition from temperature. Total sulphur and TRIS do not close the model's sulfur budget because the current model omits pyrite and other sulfur-bearing phases. Mineral abundance cannot uniquely recover porewater concentration without sample geometry, density, phase volumes, and transport history. Published mineral and isotope values therefore support comparisons and constraints, not unique estimates of unmeasured rates or mechanical parameters.

## Current inverse-use decision

No published value is currently admitted as an inverse-fit target. The two operators support a no-fit sensitivity screen only. Site-registered sample support, complete phase/species closure, and an independently justified observation-error model remain absent. Proposed mechanisms remain scenario switches and are not recast as observations. Synthetic porosity, stress, strain, permeability and rates are software benchmark values only; the source study does not provide the site-specific mechanical constraints needed to identify them.
