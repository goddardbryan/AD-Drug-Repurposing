"""Analysis-only network workflows. Python 3 reimplementation; see README limitations."""
from pathlib import Path
import argparse, json, random
from collections import defaultdict, deque
import numpy as np
import pandas as pd
import networkx as nx
ROOT=Path(__file__).resolve().parents[1]

def distances(g, targets):
    dist=dict.fromkeys(targets,0); q=deque(sorted(targets))
    while q:
        u=q.popleft()
        for v in g[u]:
            if v not in dist:dist[v]=dist[u]+1;q.append(v)
    return dist

def closest(g,a,b):
    d=distances(g,b)
    return float(np.mean([d[x] for x in a]))

def within(g,a):
    if len(a)<2:return float('nan')
    values=[]
    for u in sorted(a):
        seen={u};q=deque([(u,0)])
        while q:
            v,d=q.popleft()
            if v!=u and v in a:
                values.append(d);break
            for w in g[v]:
                if w not in seen:seen.add(w);q.append((w,d+1))
    return float(np.mean(values))


def separation(g,a,b):
    dab=(len(a)*closest(g,a,b)+len(b)*closest(g,b,a))/(len(a)+len(b))
    return dab-(within(g,a)+within(g,b))/2

def degree_pools(g,minimum):
    groups=defaultdict(list)
    for n,d in g.degree():groups[d].append(n)
    bins=[];pending=[]
    for d in sorted(groups):
        pending+=sorted(groups[d])
        if len(pending)>=minimum:bins.append(pending);pending=[]
    if pending:
        if bins:bins[-1]+=pending
        else:bins=[pending]
    return {n:i for i,b in enumerate(bins) for n in b},bins

def sample_matched(a,lookup,bins,rng):
    counts=defaultdict(int)
    for n in a:counts[lookup[n]]+=1
    return set(n for i,k in sorted(counts.items()) for n in rng.sample(bins[i],k))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('mode',choices=['separation','proximity','cluster','localization','overlap','ror'])
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--drug',help='Optional DrugBank ID; omit for all candidates')
    ap.add_argument('--permutations',type=int,default=1000)
    ap.add_argument('--seed',type=int,default=452456)
    ap.add_argument('--min-bin-size',type=int,default=2)
    args=ap.parse_args(); root=args.root; data=root/'data';out=root/'output'/args.mode;out.mkdir(parents=True,exist_ok=True)
    if args.permutations<2:ap.error('At least 2 permutations required')
    if args.mode=='ror':
        d=pd.read_csv(data/'drug_adverse_event_counts.csv'); a,b,c,e=[pd.to_numeric(d[x]) for x in 'abcd']
        valid=(a>0)&(b>0)&(c>0)&(e>0)
        d['ROR']=(a*e/(b*c)).where(valid);se=np.sqrt(1/a+1/b+1/c+1/e)
        d['CI_lower']=np.exp(np.log(d.ROR)-1.96*se);d['CI_upper']=np.exp(np.log(d.ROR)+1.96*se)
        d['inverse_signal']=(d.CI_upper<1)&(a>=5);d.to_csv(out/'results.csv',index=False);return
    edges=pd.read_csv(root/'data/ppi_interactions.tsv',sep='\t').iloc[:,[0,2]]
    g=nx.from_pandas_edgelist(edges,edges.columns[0],edges.columns[1]);g.remove_edges_from(nx.selfloop_edges(g));g=g.subgraph(max(nx.connected_components(g),key=len)).copy()
    ad=set(pd.read_csv(data/'ad_risk_genes.csv').gene.astype(str))
    modules={'AD':ad};modules.update({str(r.id):set(str(r.genes).split()) for r in pd.read_csv(data/'ad_gene_modules.csv').itertuples()})
    candidates=pd.read_csv(data/'candidate_drug_targets.csv').fillna('');controls=pd.read_csv(data/'control_drug_targets.csv').fillna('')
    if args.drug:candidates=candidates[candidates.id==args.drug]
    if candidates.empty:raise ValueError('No matching drug ID')
    universe=set(g);coverage=[]
    def mapped(name,genes):
        genes=set(genes); coverage.append({'set':name,'input_n':len(genes),'mapped_n':len(genes&universe),'missing':';'.join(sorted(genes-universe))});return genes&universe
    modules={k:mapped(k,v) for k,v in modules.items()}
    drugs={str(r.id):mapped(str(r.id),str(r.targets).split()) for r in candidates.itertuples()}
    refs={str(r.id):mapped(str(r.id),str(r.targets).split()) for r in controls.itertuples()}
    lookup,bins=degree_pools(g,args.min_bin_size);rng=random.Random(args.seed);rows=[]
    if args.mode=='cluster':
        import community
        h=g.subgraph(modules['AD']).copy()
        partition=community.best_partition(h,random_state=args.seed)
        rows=[{'gene':k,'community':v} for k,v in sorted(partition.items())]
    elif args.mode=='localization':
        for name,a in modules.items():
            if not a:continue
            observed=max(map(len,nx.connected_components(g.subgraph(a))))
            null=[max(map(len,nx.connected_components(g.subgraph(rng.sample(sorted(g),len(a)))))) for _ in range(args.permutations)]
            rows.append({'module':name,'lcc_size':observed,'null_mean':np.mean(null),'P_upper_add_one':(1+sum(v>=observed for v in null))/(args.permutations+1)})
    else:
        for name,a in drugs.items():
            for target,b in (refs if args.mode=='separation' else modules).items():
                if not a or not b:rows.append({'drug':name,'target':target,'status':'empty_mapped_set'});continue
                row={'drug':name,'target':target,'n_drug':len(a),'n_target':len(b),'status':'ok'}
                if args.mode=='separation':row['S_AB']=separation(g,a,b);row['status']='singleton_undefined' if min(len(a),len(b))<2 else 'ok'
                elif args.mode=='overlap':
                    d=distances(g,a);row.update({'overlap_n':len(a&b),'genes_distance_le_1':';'.join(sorted(x for x in b if d[x]<=1))})
                else:
                    observed=closest(g,a,b)
                    null=[closest(g,sample_matched(a,lookup,bins,rng),sample_matched(b,lookup,bins,rng)) for _ in range(args.permutations)]
                    sd=float(np.std(null));row.update({'d_c':observed,'null_mean':np.mean(null),'null_sd':sd,'Z_dc':(observed-np.mean(null))/sd if sd else float('nan'),'P_lower_add_one':(1+sum(v<=observed for v in null))/(args.permutations+1)})
                rows.append(row)
    pd.DataFrame(rows).to_csv(out/'results.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'coverage.csv',index=False)
    (out/'parameters.json').write_text(json.dumps({**vars(args),'root':str(root),'network_nodes':len(g),'network_edges':g.number_of_edges()},indent=2),encoding='utf-8')
if __name__=='__main__':main()
