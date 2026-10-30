#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict


def read_ppi(path):
    ppi = pd.read_csv(path, sep='\t').iloc[:, [0,2]].copy()
    ppi.columns = ['source','target']
    ppi = ppi.dropna()
    ppi['source'] = ppi['source'].astype(str)
    ppi['target'] = ppi['target'].astype(str)
    ppi = ppi[ppi['source'] != ppi['target']]
    ppi['a'] = ppi[['source','target']].min(axis=1)
    ppi['b'] = ppi[['source','target']].max(axis=1)
    ppi = ppi.drop_duplicates(['a','b'])[['a','b']].rename(columns={'a':'source','b':'target'})
    return ppi


def summarize_connectors(fullG, seeds):
    seeds = sorted(set(seeds) & set(fullG.nodes()))
    support = defaultdict(int)
    seed_neighbors = defaultdict(set)
    direct_edges = set()
    for i, s in enumerate(seeds):
        ns = set(fullG.neighbors(s))
        for t in seeds[i+1:]:
            if fullG.has_edge(s,t):
                direct_edges.add(tuple(sorted((s,t))))
            mids = (ns & set(fullG.neighbors(t))) - {s,t}
            for m in mids:
                support[m] += 1
                seed_neighbors[m].update([s,t])
    rows = [{'GeneSymbol':m, 'support':support[m], 'seed_neighbors':len(seed_neighbors[m])} for m in support]
    cdf = pd.DataFrame(rows).sort_values(['support','seed_neighbors','GeneSymbol'], ascending=[False,False,True]) if rows else pd.DataFrame(columns=['GeneSymbol','support','seed_neighbors'])
    return cdf, direct_edges


def build_connector_subgraph(fullG, seeds, direction_genes, top_n=30, min_seed_neighbors=8, min_support=6):
    cdf, direct_edges = summarize_connectors(fullG, seeds)
    if len(cdf):
        keep_conn = cdf[(cdf['seed_neighbors']>=min_seed_neighbors) & (cdf['support']>=min_support)].head(top_n).copy()
    else:
        keep_conn = cdf.copy()
    connectors = set(keep_conn['GeneSymbol'])
    seeds_touch = set()
    for m in connectors:
        seeds_touch.update(set(fullG.neighbors(m)) & set(seeds))
    seeds_direct = set([u for e in direct_edges for u in e])
    # Only keep nodes participating in the real connector network; do not force-add all direction genes.
    vis_nodes = connectors | seeds_touch | seeds_direct
    H = fullG.subgraph(vis_nodes).copy()
    # Remove isolates so every shown node participates in at least one real PPI edge.
    H.remove_nodes_from(list(nx.isolates(H)))
    # Keep only the largest connected component for a compact, publication-style display.
    if len(H) > 0:
        largest_component = max(nx.connected_components(H), key=len)
        H = H.subgraph(largest_component).copy()
    keep_conn = keep_conn[keep_conn['GeneSymbol'].isin(H.nodes())].copy().set_index('GeneSymbol') if len(keep_conn) else keep_conn
    return H, keep_conn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()
    root = Path(args.root).resolve()
    cfg = json.loads((root/args.config).read_text(encoding='utf-8'))
    data = root/'data'; out = root/'output/network'; out.mkdir(parents=True, exist_ok=True)
    ppi = read_ppi(data/'ppi_interactions.tsv')
    fullG = nx.Graph(); fullG.add_edges_from(ppi[['source','target']].itertuples(index=False, name=None))
    ad = pd.read_csv(data/'ad_expression_signature.csv')
    seeds = set(ad['GeneSymbol'].astype(str)) & set(fullG.nodes())
    H, keep_conn = build_connector_subgraph(fullG, seeds, seeds,
        top_n=int(cfg.get('connector_top_n',30)),
        min_seed_neighbors=int(cfg.get('connector_min_seed_neighbors',8)),
        min_support=int(cfg.get('connector_min_support',6))
    )
    nodes = pd.DataFrame({'GeneSymbol': sorted(H.nodes())})
    nodes['node_role'] = np.where(nodes['GeneSymbol'].isin(seeds), 'seed', 'connector')
    nodes['connector_support'] = nodes['GeneSymbol'].map(keep_conn['support'] if len(keep_conn) else {}).fillna(0).astype(int)
    nodes['seed_neighbor_count'] = nodes['GeneSymbol'].map(keep_conn['seed_neighbors'] if len(keep_conn) else {}).fillna(0).astype(int)
    deg = dict(H.degree())
    nodes['degree'] = nodes['GeneSymbol'].map(deg).fillna(0).astype(int)
    nodes = nodes.merge(ad[['GeneSymbol','AD_sigma','AD_direction']], on='GeneSymbol', how='left')
    nodes.loc[nodes['node_role']=='connector','AD_direction'] = 'Connector'
    nodes['AD_direction'] = nodes['AD_direction'].fillna('Not regulated / direction unavailable')
    edges = pd.DataFrame(sorted([tuple(sorted(e)) for e in H.edges()]), columns=['source','target'])
    nodes.to_csv(out/'AD_overall_connector_nodes.csv', index=False, encoding='utf-8-sig')
    edges.to_csv(out/'AD_overall_connector_edges.csv', index=False, encoding='utf-8-sig')
    comps = sorted((len(c) for c in nx.connected_components(H)), reverse=True)
    qc = pd.DataFrame([{
        'Seed_AD_signature_genes_in_PPI': len(seeds),
        'Displayed_nodes_total': len(H.nodes()),
        'Displayed_seed_nodes': int((nodes['node_role']=='seed').sum()),
        'Displayed_connector_nodes': int((nodes['node_role']=='connector').sum()),
        'Displayed_edges': len(edges),
        'Connected_components': len(comps),
        'Largest_component_nodes': comps[0] if comps else 0,
        'Connector_top_n': int(cfg.get('connector_top_n',30)),
        'Connector_min_seed_neighbors': int(cfg.get('connector_min_seed_neighbors',8)),
        'Connector_min_support': int(cfg.get('connector_min_support',6))
    }])
    qc.to_csv(out/'AD_overall_connector_QC_summary.csv', index=False, encoding='utf-8-sig')
    print(qc.to_string(index=False))

if __name__ == '__main__':
    main()
