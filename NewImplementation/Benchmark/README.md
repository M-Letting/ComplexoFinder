# Benchmark suite: thesis pipeline against new pipeline

Runs the thesis pipeline (`ComplexoFinder/`, "old") and the new one (`NewImplementation/`, "new") through one harness, on the same data, scored with the same functions.

## Running it

From the repository root:

```bash
Rscript NewImplementation/Benchmark/RunBenchmarks.R
```

or one script at a time:

| script | what it compares | data |
|---|---|---|
| [01-Discovery.R](01-Discovery.R) | discovery of discordant peptides | SWATH-MS, ProteoMaker |
| [02-Clustering.R](02-Clustering.R) | number of clusters, constraints, clustering and quantification, end to end; also the new pipeline's switches | ProteoMaker |
| [02b-CannotLinks.R](02b-CannotLinks.R) | the cannot-link rules against the truth, pair by pair | ProteoMaker |
| [03-RealData.R](03-RealData.R) | what each pipeline does on real complexes | H9 |
| [04-Report.R](04-Report.R) | builds a report from whatever results exist | |

Results are written to `NewImplementation/Benchmark/results/`, with `REPORT.md` as the place to start. `02-Clustering.R` saves its clusterings and scores them from there, so a changed metric does not need the pipelines to run again.

[BenchmarkData.R](BenchmarkData.R) holds the data loaders and [BenchmarkFunctions.R](BenchmarkFunctions.R) the wrappers around the two pipelines and the metrics.

## Data

No data is included in the repository.

- **ProteoMaker.** Simulated proteoform data: 10 conditions with 3 replicates, at three shares of regulated proteoforms (`Low`, `Med`, `High`), each with about half of the values missing (`NA`) and without missing values (`NoNA`).
  - [`Benchmark/02-DCF/00-MakeProteoMakerData.R`](../../Benchmark/02-DCF/00-MakeProteoMakerData.R) runs the six simulations with [ProteoMaker](https://github.com/computproteomics/ProteoMaker) and writes one table per dataset, with the simulation's metadata.
  - [`Benchmark/02-DCF/01-DataProcessing.R`](../../Benchmark/02-DCF/01-DataProcessing.R) prepares them: one row per peptide, a position column, median-centered samples.
  - The benchmark reads `Benchmark/02-DCF/00-Data/01-Unprepared/ProteoMaker_<dataset>.csv` and `Benchmark/02-DCF/00-Data/02-Prepared/ProteoMaker<dataset>_wide_processed.csv`.
- **SWATH-MS.** Real DIA data in which some peptides of some proteins were reduced on one day: one peptide (`1pep`), two (`2pep`), a random share (`random`) or half (`050pep`). The benchmark reads the prepared tables from `Benchmark/01-FindDCF/data/prepared/` and, as a reference, the saved results of PeCorA, ProteoForge and COPF from `Benchmark/01-FindDCF/data/results/`. The data source and the scripts that prepare it and run those tools are in the [thesis repository](https://github.com/M-Letting/MasterThesis).
- **H9.** `03-RealData.R` reads the new pipeline's saved results, so run `NewImplementation/RunComplexoFinder.R` first, and runs the thesis pipeline on its own data in `ComplexoFinder/00-Data/`.

### The truth in ProteoMaker

- A peptide's *expected profile* is the mean of the regulation profiles of the proteoforms that carry it.
- A peptide is *discordant* when its expected profile differs from the most common one in its protein.
- Peptides of a protein with the same expected profile form a *truth group*.

## 1. Discovery

- **old**: `discover_complexoforms()` of `ComplexoFinder/`, once per protein, with the pipeline's settings.
- **new**: `discover_complexoforms()` of `NewImplementation/`, one call per dataset.
- **new_global_bh**: the same with `correction = "global_bh"`.
- For SWATH-MS the saved results of PeCorA, ProteoForge and COPF are scored the same way. They are not run again.

A protein is the assembly. Every method is judged on the same peptides: all peptides of a protein with at least two. A peptide a method does not score counts as not called.

| metric | meaning |
|---|---|
| coverage | share of the peptides the method scores |
| AUROC | ranking of discordant above concordant peptides, over all peptides |
| within | the same inside each protein, pooled: can the method tell a discordant peptide from the others of its protein |
| false alarms | share called among the peptides of proteins with no discordant peptide |
| power | share of the discordant peptides called |
| false calls | share of the called peptides that are not discordant |

## 2. Clustering and 3. Quantification (ProteoMaker)

The 25 proteins with most peptides are clustered in every dataset. Each pipeline runs end to end, so the number of clusters comes from its own discovery.

| configuration | what runs |
|---|---|
| old, unconstrained | old discovery, number of clusters, VSClust |
| old (pipeline) | the same with Pearson cannot-links through `vsclust_to_restrictions()` |
| new, unconstrained | new discovery, number of clusters, VSClust |
| new (pipeline) | the same with Pearson cannot-links; the constrained run starts from the first run's centers |
| new, calibrated links | every overlapping pair tested on its own, without a tolerance (`grouping = FALSE`, `tolerance = 0`) |
| new, ZRMSD links, new, CCC links | the rules `method = "zrmsd"` and `method = "ccc"` |
| new, tolerance links (0.25) | pair by pair with a tolerance of 0.25 log2 |
| new, group links | pattern groups per position (`grouping = TRUE`), the rule the pipeline uses by default, here without a tolerance |
| new, group links (0.25) | the same with a tolerance of 0.25 log2 |
| new, scaling = standardize | VSClust on standardized instead of centered profiles |
| new, sds = global | VSClust's standard deviations from the discovery's noise model |
| new, tau-aware core | discovery with `variance_aware_core = TRUE` |
| new, global BH | discovery with `correction = "global_bh"` |
| new, true k (reference) | the new pipeline given the true number of groups |

Every clustering is scored on all peptides of the protein. A peptide is assigned to the cluster of its highest membership when that membership is at least 0.5.

| metric | meaning |
|---|---|
| k_error | difference between the clusters used and the truth groups of at least two peptides |
| assigned | share of the protein's peptides assigned to a cluster |
| pair_precision | of the peptide pairs put in one cluster, the share that belong to one truth group |
| pair_recall | of the peptide pairs in one truth group, the share put in one cluster |
| ARI | adjusted Rand index against the truth groups; each unassigned peptide is a group of its own |
| profile_R2 | share of the differences between the peptides' expected profiles that the clusters explain |
| top_fraction | of the peptides that carry proteoform 1, the largest share in one cluster |
| cl_violated | share of cannot-link pairs that sit in one cluster |
| potential | *total potential proteoforms*: a cluster's members are merged into stretches of overlapping peptides, and the cluster stands for the product over its stretches of the number of distinct proteoform IDs in each; summed over clusters. Lower is better |
| position | the same with the number of member peptides in each stretch, so it needs no truth |
| sum_unique | distinct proteoform IDs among a cluster's members, summed over clusters |
| ..._1k | the same with a flexible cutoff: a peptide is a member of every cluster it has a membership of at least 1/k in |
| ..._top | the same with every clustered peptide counted once, in the cluster of its highest membership |
| pattern, pattern_tol | potential forms with the distinct truth groups in a stretch as its alternatives; `_tol` counts expected profiles within 0.25 log2 of one another as one |
| forms present, found, recovered | truth groups of at least two peptides; clusters with at least two members; truth groups that one cluster recovers (it holds more than half of the group's peptides and more than half of its members come from the group) |
| quant_* | error (log2) of a cluster's quantified profile against the true profile of its members, both as a deviation from the protein's mean profile |

The same clustering is quantified four ways, to separate the two changes of the new pipeline: every peptide weighted by its membership (`all`) or members only (`members`), and a plain mean (`plain`) or each peptide centered on its own level first (`centred`). `quant_all_plain` is the thesis pipeline's setting and `quant_members_centred` the new one's.

## 2b. The cannot-link rules, pair by pair (ProteoMaker)

Every rule is given the same candidate pairs of overlapping peptides, and a pair truly differs when its two expected profiles do.

| metric | meaning |
|---|---|
| candidates, different | pairs a rule can link, and the share of them that truly differ |
| linked | share the rule links |
| precision | of the linked pairs, the share that truly differ |
| recall | of the truly different pairs, the share linked |
| false_links | of the pairs that do not differ, the share linked |
| AUROC | how well the rule's score ranks truly different pairs above the others |
| precision_at_n | precision of each rule's n most confident pairs, n being the number of links of the rule that makes fewest |

## 4. Real data

The H9 complexes have no ground truth, so this part describes what each pipeline does on them: peptides clustered, share called discordant, number of clusters, cannot-link pairs and the share of them in one cluster before and after the constrained run, and the share of peptides that are members.

## Things to keep in mind

- **Settings are scaled to ProteoMaker's 10 conditions.** Peptides seen in at least 4 conditions are clustered, and cannot-links need 3 shared conditions.
- **The thesis pipeline's VSClust runs are not seeded**, so its numbers can move slightly between runs.
- **The clustering, constraint and quantification steps have no real-data ground truth.** For those, the evidence is the simulation.
