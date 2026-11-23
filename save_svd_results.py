#!/usr/bin/env python3
"""Persist the refined-grid S matrices, SVD directions and null modes."""
import csv,json
from pathlib import Path
from identifiability_screen import analyze,PARAMS
OUT=Path(__file__).resolve().parent/'results'
results=[]
for law in ('kozeny_carman','power_law'):
    r=analyze(law,40,.05)
    results.append(r)
    (OUT/'identifiability_svd_refined.json').write_text(json.dumps({
      'fit_performed':False,'protocol':'identifiability_protocol.md',
      'parameter_names':[p['name'] for p in PARAMS], 'results':results},indent=2)+'\n',encoding='utf-8')
    print('persisted',r['permeability_law'],'rank',r['matrix_rank'],'distinguishable',r['distinguishable_direction_count'],flush=True)
# 20-cell sensitivity matrices were run separately as the coarse-grid numerical check.
coarse_values=[3.368122033380916,0.06614410505900903,3.309055245977657e-05,0.0]
rows=[]
for law in ('Kozeny-Carman','Power law'):
    rows.append({'permeability_law':law,'nx':20,'dt_s':0.2,'matrix_rank_rel_1e-8':3,
                 'distinguishable_directions_sigma_ge_1':1,'singular_values':';'.join(f'{v:.10g}' for v in coarse_values)})
for r in results:
    rows.append({'permeability_law':r['permeability_law'],'nx':40,'dt_s':0.05,'matrix_rank_rel_1e-8':r['matrix_rank'],
                 'distinguishable_directions_sigma_ge_1':r['distinguishable_direction_count'],
                 'singular_values':';'.join(f'{v:.10g}' for v in r['singular_values'])})
with (OUT/'identifiability_rank_summary.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
