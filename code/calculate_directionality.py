"""Recompute supported directionality statistics; preserve provenance of stored Zdc."""
from pathlib import Path
import csv
import math
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/directionality'
OUT.mkdir(parents=True, exist_ok=True)

def read(name):
    with (ROOT / 'data' / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

ad = read('ad_expression_signature.csv')
lincs = read('lincs_drug_gene_zscores.csv')
assert len({r['GeneSymbol'] for r in ad}) == len(ad) == 291
assert len({r['GeneSymbol'] for r in lincs}) == len(lincs)
sigma = {r['GeneSymbol']: float(r['AD_sigma']) for r in ad}
drugs = list(lincs[0])[1:]
assert len(drugs) == 8
stats, evidence = {}, []
for drug in drugs:
    values = []
    for row in lincs:
        gene = row['GeneSymbol']
        if gene not in sigma or not row[drug]:
            continue
        value = float(row[drug])
        if not math.isfinite(value):
            continue
        product = sigma[gene] * value
        relation = 'reversed' if product < 0 else 'reinforced' if product > 0 else 'unchanged'
        values.append(relation)
        evidence.append(dict(Drug=drug, GeneSymbol=gene, AD_sigma=sigma[gene], drug_z=value, relation=relation))
    n = len(values)
    rev, rein = values.count('reversed'), values.count('reinforced')
    pad = (rev - rein) / n
    z = pad * math.sqrt(n)
    stats[drug] = dict(N_used=n, N_reversed=rev, N_reinforced=rein, N_unchanged=n-rev-rein, pAD=pad, pAD_Z=z)


for filename, records in [('overall_directionality.csv', [dict(Drug=k, **v) for k,v in stats.items()]), ('gene_level_directionality.csv', evidence)]:
    with (OUT/filename).open('w', encoding='utf-8-sig', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
print('Directionality calculation completed.')
