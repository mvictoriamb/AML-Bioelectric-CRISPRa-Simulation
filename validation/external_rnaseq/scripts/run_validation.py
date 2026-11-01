"""Reproduce association replication from downloaded cBioPortal data, not Keras predictions.
Run: python run_validation.py --data /path/to/external_rnaseq --discovery /path/to/Leukemia_GSE9476.csv --output results
Requires numpy, scipy, matplotlib. Full cohort matrices remain outside deliverables.
"""
import argparse,csv,hashlib,json,zipfile
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import spearmanr,pearsonr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def clinical(path):
    with path.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))

def expression(path):
    found={}
    with path.open(encoding='utf-8-sig') as f:
        r=csv.reader(f,delimiter='\t'); header=next(r)
        for row in r:
            if row[0] in ('CREB5','KCNJ15','LOC401317'):
                assert row[0] not in found
                found[row[0]]=dict(zip(header[2:],map(float,row[2:])))
    assert 'CREB5' in found and 'KCNJ15' in found
    return found

def bh(p):
    p=np.array(p); order=np.argsort(p); n=len(p); q=np.empty(n)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]); return q

def rho_boot(x,y,rng,b=10000):
    out=[]; n=len(x)
    for start in range(0,b,500):
        ix=rng.integers(0,n,(min(500,b-start),n)); a=rankdata(x[ix],axis=1); c=rankdata(y[ix],axis=1)
        a-=a.mean(axis=1,keepdims=True); c-=c.mean(axis=1,keepdims=True)
        out.extend(np.sum(a*c,axis=1)/np.sqrt(np.sum(a*a,axis=1)*np.sum(c*c,axis=1)))
    return np.array(out)

def perm_p(x,y,rng,b=99999):
    a=rankdata(x);a-=a.mean();c=rankdata(y);c-=c.mean()
    denom=np.linalg.norm(a)*np.linalg.norm(c); obs=np.dot(a,c)/denom; extreme=0
    for start in range(0,b,500):
        perm=rng.permuted(np.tile(c,(min(500,b-start),1)),axis=1)
        extreme+=np.count_nonzero(np.abs(perm@a/denom)>=abs(obs)-1e-14)
    return float((extreme+1)/(b+1))

def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main(a):
    a.output.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(42)
    cohort_rows={}; flow={}; inputs=[]
    # Clinical rules fixed before correlation calculation. No selection by expression.
    for name,d,filename in [('TCGA-LAML',a.data,'data_mrna_seq_v2_rsem.txt'),
                            ('Beat-AML',a.data/'beataml','data_mrna_seq_rpkm.txt')]:
        ex=expression(d/filename); sm=clinical(d/'data_clinical_sample.txt'); patients=clinical(d/'data_clinical_patient.txt')
        pat={p['PATIENT_ID']:p for p in patients}
        candidates=[]
        for s in sm:
            sid=s['SAMPLE_ID']; pid=s['PATIENT_ID']
            if sid not in ex['CREB5'] or sid not in ex['KCNJ15']: continue
            if name=='TCGA-LAML': eligible=s['SAMPLE_TYPE']=='Primary'
            else: eligible=s['SAMPLE_TIMEPOINT']=='Denovo' and pat[pid]['DIAGNOSIS']=='Acute Myeloid Leukaemia (AML) and Related Precursor Neoplasms'
            if not eligible: continue
            x,y=ex['CREB5'][sid],ex['KCNJ15'][sid]
            if not np.isfinite(x) or not np.isfinite(y): continue
            candidates.append({**s,'CREB5':x,'KCNJ15':y})
        bypatient=defaultdict(list)
        for s in candidates: bypatient[s['PATIENT_ID']].append(s)
        chosen=[]
        for pid,ss in sorted(bypatient.items()):
            # Prefer marrow then stable sample ID among eligible same-timepoint samples.
            ss.sort(key=lambda s:(0 if s.get('SAMPLE_SITE','')=='Bone Marrow Aspirate' else 1,s['SAMPLE_ID']))
            s=ss[0]; p=pat[pid]
            chosen.append({'cohort':name,'patient_id':pid,'sample_id':s['SAMPLE_ID'],
                'CREB5':s['CREB5'],'KCNJ15':s['KCNJ15'],
                'site':s.get('SAMPLE_SITE',s.get('PRIMARY_SITE','')),
                'age':p.get('AGE',p.get('AGE_AT_DIAGNOSIS','')),
                'sex':p.get('SEX',''),'FAB':s.get('FAB',p.get('FAB',''))})
        assert len({s['patient_id'] for s in chosen})==len(chosen)
        cohort_rows[name]=chosen
        flow[name]={'clinical_samples':len(sm),'expression_columns':len(ex['CREB5']),
            'eligible_complete_samples':len(candidates),'unique_patients':len(chosen),
            'patients_with_multiple_eligible_samples':sum(len(v)>1 for v in bypatient.values()),
            'LOC401317_row_present':'LOC401317' in ex,'sites':dict(Counter(s['site'] for s in chosen))}
        for f in (filename,'data_clinical_sample.txt','data_clinical_patient.txt','meta_study.txt'):
            path=d/f
            inputs.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size})
    with a.discovery.open(encoding='utf-8-sig') as f: micro=list(csv.DictReader(f))
    discovery=np.array([[float(s['205931_s_at']),float(s['210119_at'])] for s in micro if s['type']=='AML'])
    discovery_boot=rho_boot(discovery[:,0],discovery[:,1],rng)
    results=[]; boot={}; sensitivities=[]
    # Family: four external tests, two cohorts x Spearman primary / Pearson secondary.
    for name,rows in cohort_rows.items():
        xy=np.array([[r['CREB5'],r['KCNJ15']] for r in rows]); x,y=xy.T
        rho=float(spearmanr(x,y).statistic); boot[name]=rho_boot(x,y,rng)
        ci=np.quantile(boot[name],[.025,.975]); diff=np.quantile(boot[name]-discovery_boot,[.025,.975])
        results.append({'cohort':name,'test':'Spearman_primary','n':len(x),'effect':rho,
            'ci95_low':float(ci[0]),'ci95_high':float(ci[1]),'p_raw':perm_p(x,y,rng),
            'scale':'ranks of provided gene expression','difference_vs_discovery_rho':rho-float(spearmanr(*discovery.T).statistic),
            'difference_ci95_low':float(diff[0]),'difference_ci95_high':float(diff[1])})
        xx,yy=(np.log2(x+1),np.log2(y+1)) if name=='TCGA-LAML' else (x,y)
        pear=pearsonr(xx,yy); pci=pear.confidence_interval(.95)
        results.append({'cohort':name,'test':'Pearson_secondary','n':len(x),'effect':float(pear.statistic),
            'ci95_low':float(pci.low),'ci95_high':float(pci.high),'p_raw':float(pear.pvalue),
            'scale':'log2(RSEM+1)' if name=='TCGA-LAML' else 'as provided (includes negative values; not logged again)',
            'difference_vs_discovery_rho':'','difference_ci95_low':'','difference_ci95_high':''})
        for site in sorted({r['site'] for r in rows}):
            subset=[r for r in rows if r['site']==site]
            if len(subset)>=10:
                sensitivities.append({'cohort':name,'site':site,'n':len(subset),
                    'rho_descriptive':float(spearmanr([r['CREB5'] for r in subset],[r['KCNJ15'] for r in subset]).statistic)})
    for row,q in zip(results,bh([r['p_raw'] for r in results])):row['q_BH_family4']=float(q)
    write_csv(a.output/'external_tests.csv',results)
    write_csv(a.output/'selected_samples_expression.csv',[s for rows in cohort_rows.values() for s in rows])
    write_csv(a.output/'site_sensitivity_descriptive.csv',sensitivities)
    report={'mapping':'CREB5 (9586), not independently resolved CREB5-AS1; historical probe 205931_s_at',
        'flow':flow,'results':results,'discovery_AML':{'n':len(discovery),
        'rho':float(spearmanr(*discovery.T).statistic),'ci95_bootstrap':np.quantile(discovery_boot,[.025,.975]).tolist()},
        'seed':42,'bootstrap_replicates':10000,'permutations':99999,
        'BH_family':'4 external tests: 2 cohorts x 2 methods, no discovery tests included',
        'datahub_commit':(a.data/'datahub_commit.txt').read_text().strip(),
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__},'inputs':inputs,
        'predictive_model_evaluation':'not performed: missing serialized preprocessing and cross-platform transport not defined'}
    (a.output/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
    plot=[('GSE9476: sólo AML',discovery,float(spearmanr(*discovery.T).statistic),'microarray, escala original')]
    for name,rows in cohort_rows.items():
        xy=np.array([[r['CREB5'],r['KCNJ15']] for r in rows])
        if name=='TCGA-LAML': xy=np.log2(xy+1)
        rho=next(r['effect'] for r in results if r['cohort']==name and r['test']=='Spearman_primary')
        plot.append((name,xy,rho,'log2(RSEM+1)' if name=='TCGA-LAML' else 'escala suministrada'))
    for ax,(title,xy,rho,scale) in zip(axes,plot):
        ax.scatter(*xy.T,s=17,alpha=.6,color='#286b8e'); ax.set_title(f'{title}\nn={len(xy)}; Spearman ρ={rho:.3f}')
        ax.set_xlabel('205931_s_at' if title.startswith('GSE') else 'CREB5');ax.set_ylabel('210119_at' if title.startswith('GSE') else 'KCNJ15');ax.grid(alpha=.2)
        ax.text(.02,.98,scale,transform=ax.transAxes,va='top',fontsize=8)
    fig.suptitle('Réplica de asociación: no equivale a validar el clasificador ni una terapia')
    fig.savefig(a.output/'replication.png',dpi=180);fig.savefig(a.output/'replication.svg');plt.close(fig)
    print(json.dumps({k:report[k] for k in ('flow','results','discovery_AML')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('data','discovery','output'): p.add_argument('--'+arg,type=Path,required=True)
    main(p.parse_args())
