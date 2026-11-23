#!/usr/bin/env python3
"""Expand saved 40-cell SVD matrices into full mode and nullspace tables."""
import csv,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parent; out=root/'results'
p=out/'identifiability_svd_refined.json'; j=json.loads(p.read_text())
params=j['parameter_names']; modes=[]; nulls=[]
for r in j['results']:
    S=np.asarray(r['S_whitened'],float); U,s,Vt=np.linalg.svd(S,full_matrices=True)
    rank=int(np.sum((s/s[0])>=1e-8)) if s[0] else 0
    r['full_right_singular_vectors']=Vt.tolist()
    r['numerical_rank_relative_cutoff_1e-8']=rank
    r['practical_rank_sigma_ge_1']=int(np.sum(s>=1.0))
    r['relative_condition_number_nonzero']=float(s[0]/s[rank-1]) if rank else None
    for i,val in enumerate(s):
        vec=Vt[i]
        modes.append({'permeability_law':r['permeability_law'],'mode':i+1,'singular_value':float(val),
                      'relative_to_largest':float(val/s[0]) if s[0] else 0.0,
                      'distinguishable_sigma_ge_1':bool(val>=1.0),
                      'right_singular_vector':'; '.join(f'{params[k]}={vec[k]:+.6f}' for k in np.argsort(np.abs(vec))[::-1])})
    for i in range(rank,len(params)):
        vec=Vt[i]
        top=np.argsort(np.abs(vec))[::-1][:5]
        old=r['null_directions'][i-rank] if i-rank<len(r['null_directions']) else {}
        nulls.append({'permeability_law':r['permeability_law'],'null_mode':i-rank+1,
                      'dominant_combination':'; '.join(f'{params[k]}={vec[k]:+.4f}' for k in top),
                      'measurement_to_break_direction':old.get('measurement_to_break_direction','See parameter-specific measurement recommendations in report.')})
(out/'identifiability_svd_refined.json').write_text(json.dumps(j,indent=2)+'\n')
for name,rows in [('identifiability_modes.csv',modes),('identifiability_nullspace.csv',nulls)]:
    with (out/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('wrote full vectors and',len(nulls),'law-specific null basis rows')
