"""Download the exact public cBioPortal matrices used in this analysis.
Python standard library only. About 156 MB plus metadata. No authentication.
Usage: python download_data.py --destination data_download
"""
import argparse,json,hashlib,urllib.request
from pathlib import Path
COMMIT='e1b75cd25bf7436d16719b7d76117ad8f3f66e54'

def main(dest):
    sources=[]
    for study,sub,expr in [('laml_tcga','', 'data_mrna_seq_v2_rsem.txt'),
                           ('aml_ohsu_2018','beataml','data_mrna_seq_rpkm.txt')]:
        folder=dest/sub;folder.mkdir(parents=True,exist_ok=True)
        for filename in (expr,'data_clinical_sample.txt','data_clinical_patient.txt','meta_study.txt',
                         'meta_mrna_seq_v2_rsem.txt' if not sub else 'meta_mrna_seq_rpkm.txt'):
            host='raw.githubusercontent.com/cBioPortal/datahub' if filename.startswith('meta_') else 'media.githubusercontent.com/media/cBioPortal/datahub'
            url=f'https://{host}/{COMMIT}/public/{study}/{filename}'
            target=folder/filename
            urllib.request.urlretrieve(url,target)
            sources.append({'url':url,'file':str(target),'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
            print(study,filename,target.stat().st_size)
    (dest/'datahub_commit.txt').write_text(COMMIT+'\n')
    (dest/'download_manifest.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--destination',type=Path,default=Path('data_download'))
    main(p.parse_args().destination)
