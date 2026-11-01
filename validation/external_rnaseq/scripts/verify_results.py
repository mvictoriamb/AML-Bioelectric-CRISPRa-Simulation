"""Offline numerical checks of the published extract and testing family.
Does not repeat download, cohort selection, bootstrap or permutation tests.
Run from any working directory: python scripts/verify_results.py
"""
import csv,json
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr,spearmanr,false_discovery_control

root=Path(__file__).resolve().parents[1]
report=json.loads((root/'results/results.json').read_text(encoding='utf-8'))
with (root/'data/selected_samples_expression.csv').open(encoding='utf-8',newline='') as f:
    samples=list(csv.DictReader(f))
tests=report['results']
assert len(tests)==4
assert np.allclose(false_discovery_control([r['p_raw'] for r in tests],method='bh'),
                   [r['q_BH_family4'] for r in tests],rtol=1e-12,atol=0)
for cohort in ('TCGA-LAML','Beat-AML'):
    rows=[r for r in samples if r['cohort']==cohort]
    assert len(rows)==len({r['patient_id'] for r in rows})
    x=np.array([float(r['CREB5']) for r in rows]); y=np.array([float(r['KCNJ15']) for r in rows])
    primary=next(r for r in tests if r['cohort']==cohort and r['test']=='Spearman_primary')
    assert primary['n']==len(rows)
    assert abs(spearmanr(x,y).statistic-primary['effect'])<1e-12
    if cohort=='TCGA-LAML': x,y=np.log2(x+1),np.log2(y+1)
    secondary=next(r for r in tests if r['cohort']==cohort and r['test']=='Pearson_secondary')
    result=pearsonr(x,y)
    assert abs(result.statistic-secondary['effect'])<1e-12
    assert np.isclose(result.pvalue,secondary['p_raw'],rtol=1e-10,atol=0)
with (root/'data/discovery_aml_pair.csv').open(encoding='utf-8',newline='') as f:
    discovery=list(csv.DictReader(f))
assert len(discovery)==26 and all(r['type']=='AML' for r in discovery)
rho=spearmanr([float(r['205931_s_at']) for r in discovery],
              [float(r['210119_at']) for r in discovery]).statistic
assert abs(rho-report['discovery_AML']['rho'])<1e-12
print('PASS: patient uniqueness, external Spearman/Pearson, discovery subset and BH family of four tests.')
