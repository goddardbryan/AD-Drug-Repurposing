"""Tanimoto from externally generated PaDEL MACCS binary fingerprints; no substitute fingerprint engine."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--fingerprints',required=True,help='CSV: id followed by MACCS binary columns, same encoding for all drugs')
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args()
    f=pd.read_csv(a.fingerprints).set_index('id')
    if not f.index.is_unique or not np.isin(f.to_numpy(),[0,1]).all():raise ValueError('Unique IDs and binary bits required')
    candidates=pd.read_csv(a.root/'data/candidate_drug_targets.csv');controls=pd.read_csv(a.root/'data/control_drug_targets.csv');rows=[]
    for c in candidates.id:
        for r in controls.id:
            if c not in f.index or r not in f.index:rows.append({'candidate':c,'control':r,'T':np.nan,'status':'missing_fingerprint'});continue
            x=f.loc[c].to_numpy().astype(bool);y=f.loc[r].to_numpy().astype(bool);n=np.sum(x|y)
            rows.append({'candidate':c,'control':r,'T':np.sum(x&y)/n if n else np.nan,'status':'ok' if n else 'empty_fingerprint'})
    out=a.root/'output/structure';out.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(out/'tanimoto.csv',index=False)
if __name__=='__main__':main()
