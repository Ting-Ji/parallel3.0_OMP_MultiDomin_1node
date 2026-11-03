from pathlib import Path
import json,csv,math
import numpy as np
R=Path(__file__).resolve().parent/'full_rod_beta06_axial'
out={}
for name in ['old3','exact']:
 p=R/'runs'/name/'DBEM1/output/rod_front_case_56/validation_state.csv'
 with p.open() as f:a=list(csv.DictReader(f))
 x=np.array([float(r['x']) for r in a]);t=np.array([int(r['step'])*.0021 for r in a]);u=np.array([float(r['ux']) for r in a]);d=(6-x)/50
 truth=np.maximum(t-d,0)/400;refnorm=np.linalg.norm(truth)
 shifts=np.linspace(-2*.0021,2*.0021,401)
 fits=[]
 for shift in shifts:
  ref=np.maximum(t-shift-d,0)/400;scale=float(u@ref/(ref@ref));fits.append((float(np.linalg.norm(u-scale*ref)/refnorm),float(shift),scale))
 tip=[r for r in a if float(r['x'])>5.999]
 out[name]={'final_tip':{k:tip[-1][k] for k in ['step','x','ux','uy','uz','tx','ty','tz']},'final_tip_reference':.1176/400,'best_fitted_shift_and_scale_diagnostic_only':min(fits),'applied_axial_load_audited':all(abs(float(r['tx'])-1)<1e-10 and abs(float(r['ty']))<1e-10 and abs(float(r['tz']))<1e-10 for r in tip if int(r['step'])>0)}
(R/'reference_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
