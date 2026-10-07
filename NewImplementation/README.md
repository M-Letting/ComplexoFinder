# ComplexoFinder: the new implementation

Finds complexoforms in peptide-level proteomics data: groups of peptides of a protein complex that change together across the conditions, apart from the rest of the complex.

This folder holds the reworked pipeline. [`ComplexoFinder/`](../ComplexoFinder/) holds the pipeline of the MSc thesis, which the [benchmark suite](Benchmark/README.md) compares it with.

## Running it

From the repository root:

```bash
Rscript NewImplementation/01-DataProcessing/ProcessData.R   # once; builds NewImplementation/00-Data/
Rscript NewImplementation/RunComplexoFinder.R
```

- The settings (dataset, complexes, switches) are at the top of [RunComplexoFinder.R](RunComplexoFinder.R).
- Results go to `NewImplementation/10-Results/<dataset>/`.
- The discovery step is saved and reused while the data and its settings stay the same.

## Input

No data is included in the repository.

- **Peptide tables.** [ProcessData.R](01-DataProcessing/ProcessData.R) reads one CSV per data type from `ComplexoFinder/00-Data/`:
  `AcetylatedLysines_Peptide.csv`, `Deglycosylated_Peptide.csv`, `FreeCysteines_Peptide.csv`, `NonModified_Peptide.csv`, `NonModified_Protein.csv`, `Phospho_Peptide.csv` and `ReversiblyModifiedCysteines_Peptide.csv`.
  They are made from the raw exports by [`ComplexoFinder/01-DataProcessing/DataLoadPreprocess.R`](../ComplexoFinder/01-DataProcessing/DataLoadPreprocess.R).
  The intensity columns are named `<cell line>_B<batch>_D<day>`, for example `H9_B1_D17`.
- **Complexes.** The human table of the [Complex Portal](https://www.ebi.ac.uk/complexportal/) is downloaded on first use and cached in `NewImplementation/00-LookupTables/`.
- **Sequences (optional).** The sequence plots read `ComplexoFinder/00-LookupTables/all_human_sequences.rds`, built by [`ComplexoFinder/00-LookupTables/01-GetLookupTables.R`](../ComplexoFinder/00-LookupTables/01-GetLookupTables.R). Without it they fetch the sequences from UniProt.

## What the pipeline does

1. **Lookup** ([GetLookupTables.R](00-LookupTables/GetLookupTables.R)). Builds the complex lookup from the Complex Portal's expanded participant list, so a complex nested in another is replaced by its proteins. Isoform and processed-chain identifiers are collapsed to the canonical accession, in the lookup and in the data.
2. **Data** ([ProcessData.R](01-DataProcessing/ProcessData.R)). Splits the tables by cell line, sets zeros to missing, log2-transforms and median-centers every sample. No completeness filter is applied and no value is imputed.
3. **Complexes** ([IdentifyComplexes.R](02-Identification/IdentifyComplexes.R)). Matches peptides to complexes and resolves peptides shared between proteins. A complex is kept when at least `min_found` of its members are detected and at least `min_fraction` of them; a set of alternative members (paralogs listed for one position) counts as one member, detected when any of them is. Every peptide gets an identifier that is unique within its complex.
4. **Discovery** ([DiscoverComplexoforms.R](02-Identification/DiscoverComplexoforms.R)). One call over all kept complexes finds the discordant peptides: those whose change across the conditions differs from their complex's.
   - A peptide is tested when it has a value in at least `min_conditions` conditions.
   - *Replicate noise.* A peptide's variance between replicates is pooled over the conditions and shrunk towards a prior fitted per data type (limma).
   - *Centering.* Each peptide is centered on its own mean, which leaves its profile.
   - *Coherent core and reference.* Within a complex the peptides are clustered by profile (average linkage, cut at `core_z` in units of replicate noise); the largest cluster is the core. A peptide's reference is the mean profile of the core, without the peptide itself.
   - *Between-peptide variance.* How much peptides of one form differ beyond replicate noise, estimated once for the dataset by matching the median standardized residual about the references.
   - *Test.* The squared deviations from the reference, each divided by its variance (replicate noise, the uncertainty of the reference and the between-peptide variance), are summed and compared with a chi-square distribution.
   - *Multiple testing.* Bonferroni within the complex, then Benjamini-Hochberg across all complexes.
   - *Groups.* The discordant peptides of a complex are clustered on their profile distance in units of noise (Ward linkage, dynamic tree cut with a minimum cluster size of one). A group of at least two peptides is a differential complexoform (`dCF1`, `dCF2`, ...); a discordant peptide on its own gets `dCF-1`. The number of clusters for the next steps is one plus the number of groups.
5. **Cannot-link constraints** ([Constraints.R](03-Constraints/Constraints.R)). Overlapping peptides of a protein are versions of one stretch of it. Versions that change alike may share a complexoform; versions that differ may not. The overlapping peptides of a stretch are grouped by profile: they are joined into a tree, and a branching splits when the profiles of its two branches differ by more than noise and a tolerance (tested against the noise model of the discovery, Benjamini-Hochberg over all branchings). Peptides in different groups get a cannot-link.
6. **Clustering** ([VsClust.R](04-Clustering/VsClust.R)). VSClust (variance-sensitive fuzzy clustering) on the centered profiles of the peptides seen in at least `min_conditions_clustering` conditions. It runs twice: without constraints, and then, started from the centers of the first run, with the restrictions that the cannot-links and the first run give (a linked peptide is kept out of its partner's cluster).
7. **Quantification** ([Quantification.R](05-Quantification/Quantification.R)). A complexoform's abundance per sample is the membership-weighted mean of its members (membership of at least `membership_cutoff`), each centered on its own level first, with the group's level added back.
8. **Plots** ([09-Visualization/](09-Visualization/)). Profiles and memberships, quantified complexoforms, peptides along the sequence, and a PCA of the peptides.

## Settings

All are set at the top of [RunComplexoFinder.R](RunComplexoFinder.R).

| setting | default | meaning |
|---|---|---|
| `dataset` | `"H9"` | cell line to analyse (`"H9"` or `"IMR90"`) |
| `complexes` | a list of names | the complexes that are clustered and quantified; discovery always runs on all kept complexes |
| `include_NM` | `FALSE` | include the protein-level table |
| `min_found`, `min_fraction` | 2, 0.75 | detected members a complex needs, as a number and as a share |
| `alpha` | 0.05 | threshold on the adjusted p-value of the discovery |
| `min_conditions`, `min_reps_per_condition` | 2, 1 | conditions a peptide needs to be tested, and replicates a condition needs to count |
| `core_z` | 1.5 | cut height of the coherent core |
| `tau_aware_core` | `FALSE` | `TRUE`: find the core a second time with the between-peptide variance in its distance |
| `correction` | `"two_stage"` | `"global_bh"`: Benjamini-Hochberg across all peptides only |
| `deep_split` | 2 | 0 to 4: sensitivity of the dynamic tree cut that groups the discordant peptides |
| `reuse_discovery` | `TRUE` | reuse the saved discovery when data and settings match |
| `min_conditions_clustering` | 10 | conditions a peptide needs to be clustered |
| `scaling` | `"center"` | `"standardize"`: also scale every profile to unit variance, so only its shape counts |
| `sds_source` | `"complex"` | where VSClust's standard deviations come from: limma within the complex, or `"global"` for the discovery's noise model |
| `seed` | 1 | seed of VSClust's random starts |
| `include_constraints` | `TRUE` | use cannot-link constraints |
| `cannot_link_grouping` | `TRUE` | `FALSE`: test every overlapping pair on its own instead of grouping per position |
| `cannot_link_tolerance` | `NULL` | difference between two profiles (log2) that does not count; `NULL` uses the between-peptide standard deviation of the discovery, 0 links on any difference beyond noise |
| `cannot_link_th` | `NULL` | adjusted p-value for a link; `NULL` is 0.05 |
| `min_shared_cond` | 10 | conditions two overlapping peptides must share to be compared |
| `membership_cutoff` | 0.5 | membership a peptide needs to be a member of a complexoform |

`create_cannot_link()` also has the rules of the thesis pipeline (`method = "pearson"`, `"zrmsd"`, `"ccc"`), which the benchmark uses.

## Output

In `NewImplementation/10-Results/<dataset>/`:

- `complex_pruning_summary.csv`: the complexes and whether they were kept.
- `discovery.rds`, `discovery_peptides.csv`, `discovery_summary.csv`: the discovery over all kept complexes.
- `run_summary.csv`: one row per analysed complex, including whether the cannot-links are kept and whether VSClust converged (it stops after 1000 iterations and the pipeline warns when that limit is reached).
- One folder per analysed complex with `complexoform_results.rds` (everything from every step), `cluster_summary.csv`, `complexoform_abundance.csv`, `complexoform_group_stats.csv` and the plots.

## Requirements

R with `data.table`, `here`, `limma`, `vsclust`, `dynamicTreeCut` and `matrixStats`; for the plots `ggplot2`, `cowplot`, `scico`, `scales`, `ComplexHeatmap`, `circlize` and `pcaMethods`.

## Benchmark

[Benchmark/](Benchmark/README.md) runs this pipeline and the thesis pipeline through one harness, on data with a known truth.
