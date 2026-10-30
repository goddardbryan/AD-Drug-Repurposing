# AD Drug Repurposing Analysis

Computational workflows for Alzheimer’s disease (AD) drug prioritization using disease semantic similarity, target-network separation, network proximity, community detection, pharmacovigilance statistics, and transcriptional directionality.

## Repository structure

```text
code/               Analysis scripts
data/               PPI network and processed analysis inputs
config.json         Connector-network parameters
requirements.txt    Python dependencies
```

## Installation

Create a Python environment and install the dependencies:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
```

The Python workflows were checked with NumPy 1.21.5, pandas 1.4.4, NetworkX 2.7, SciPy 1.9.1, and python-louvain 0.16.

Disease semantic similarity additionally requires R and the Bioconductor package [meshes](https://bioconductor.org/packages/meshes/):

```r
install.packages("BiocManager")
BiocManager::install("meshes")
```

## Usage

Run commands from the repository root. Generated tables and run parameters are written to `output/`.

### Network analyses

```bash
python code/network_analysis.py separation
python code/network_analysis.py proximity --permutations 1000 --seed 452456
python code/network_analysis.py cluster
python code/network_analysis.py localization --permutations 1000 --seed 452456
python code/network_analysis.py overlap
```

| Mode | Analysis |
| --- | --- |
| `separation` | Target-network separation between candidate and control drugs |
| `proximity` | Degree-matched drug-to-module network proximity |
| `cluster` | Louvain communities in the AD-induced network |
| `localization` | Largest connected component enrichment against random gene sets |
| `overlap` | Shared targets and module genes within one network edge of drug targets |

Add `--drug DB01254` to restrict `separation`, `proximity`, or `overlap` to dasatinib. Omitting this argument processes all candidate drugs. Complete pairwise and permutation analyses may take substantial time. Use `--help` for additional arguments.

### Disease semantic similarity

```bash
Rscript code/disease_semantic_similarity.R .
```

The workflow calculates Wang similarity between AD (MeSH `D000544`) and the supplied disease identifiers using `mesh_semantic_data_ah107105.rds`. Scores are joined to drug indications, with `DSS >= 0.5` used for selection.

### Chemical structure similarity

Generate MACCS fingerprints using the PaDEL configuration corresponding to the intended analysis, then run:

```bash
python code/structure_similarity.py --fingerprints /path/to/padel_maccs.csv
```

The fingerprint CSV must contain an `id` column with unique DrugBank identifiers followed by consistently ordered binary MACCS columns. The script calculates Tanimoto similarity as the number of shared active bits divided by the number of active bits in their union. Missing or empty fingerprints are reported explicitly. SMILES are provided in the candidate and control input tables; PaDEL fingerprints and configuration are external prerequisites.

### Reporting odds ratios

```bash
python code/network_analysis.py ror
```

Calculates reporting odds ratios and log-Wald 95% confidence intervals from the supplied four-cell counts. An inverse reporting signal requires an upper confidence limit below 1 and at least five reported events. Nonpositive counts yield missing estimates.

### Connector network and transcriptional directionality

```bash
python code/build_connector_network.py
python code/calculate_directionality.py
```

The connector workflow uses the parameters in `config.json` and exports node, edge, and quality-control tables. The directionality workflow combines AD expression directions with LINCS z-scores to calculate reversal and reinforcement counts, pAD, and pAD_Z. Gene-level directionality can be joined to the connector-node table using `GeneSymbol`.

## Input data

All input files are stored together in `data/`. The supplementary-table inputs are listed below.

| File | Source and contents |
| --- | --- |
| `candidate_drug_targets.csv` | Table S1 B: candidate identifiers, targets, and SMILES |
| `control_drug_targets.csv` | Table S1 C: control identifiers, targets, and SMILES |
| `disease_mesh_ids.csv` | Table S2 A: disease names and MeSH identifiers |
| `ad_risk_genes.csv` | Table S3 A: AD-associated genes |
| `ad_gene_modules.csv` | Table S3 B: module identifiers and gene sets |
| `drug_indications.csv` | Table S4 E: drug–indication relationships |
| `drug_adverse_event_counts.csv` | Table S5 A: drug identifiers and four-cell counts |
| `mesh_semantic_data_ah107105.rds` | Annotation object associated with the original semantic-similarity workflow |

Additional inputs in `data/` are `ad_expression_signature.csv` (AD expression directions), `lincs_drug_gene_zscores.csv` (LINCS z-scores for 195 AD-signature genes across eight drugs), and `ppi_interactions.tsv` (PPI edges). The PPI table uses the first and third columns as interaction endpoints.

## Methods and reproducibility

- **Network preprocessing:** interactions are treated as undirected; self-loops are removed and the largest connected component is retained. Network workflows report mapped gene counts and missing identifiers.
- **Target separation:** cross-set distance is the pooled mean of nearest distances from A to B and B to A. S_AB subtracts half the sum of the two within-set nearest-neighbor means, excluding self-distances. Within-set distance for a singleton is undefined and produces a missing separation estimate.
- **Proximity:** the observed distance is the mean nearest distance from drug targets to module genes. Random sets preserve set sizes and degree-bin membership. Sampling is without replacement within bins and may include original nodes. Z_dc uses the population standard deviation of the null distances; zero variance produces a missing Z-score. 
- **Localization:** the observed induced-subgraph component size is compared with uniformly sampled gene sets of equal size. The empirical upper-tail probability is `(k + 1) / (N + 1)`.
- **Clustering:** Louvain uses a fixed seed on the AD-induced network. Recomputed communities are distinct from the previously defined modules in Table S3. Proximity analyses use those supplied modules and retain their A–E identifiers.
- **Pharmacovigilance:** ROR calculations use the supplied aggregate four-cell counts.

## References and attribution

Network workflows were adapted from the supplied proximity, separation, localization, and Louvain analysis scripts. The upstream network-method implementation is available in [emreg00/toolbox](https://github.com/emreg00/toolbox). Disease semantic similarity uses [meshes](https://bioconductor.org/packages/meshes/).
