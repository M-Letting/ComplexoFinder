# ComplexoFinder

Finds complexoforms in peptide-level proteomics data: groups of peptides of a protein complex that change together across the conditions, apart from the rest of the complex.

## Contents

| folder | what it holds |
|---|---|
| [NewImplementation/](NewImplementation/README.md) | The pipeline. Start here. |
| [NewImplementation/Benchmark/](NewImplementation/Benchmark/README.md) | The benchmark suite: the pipeline against the thesis pipeline, on data with a known truth. |
| [ComplexoFinder/](ComplexoFinder/) | The pipeline of the MSc thesis, kept as the baseline the benchmark compares with. |
| [Benchmark/02-DCF/](Benchmark/02-DCF/) | The two scripts that generate and prepare the simulated ProteoMaker datasets of the benchmark. |

The repository holds code only: no data, lookup tables, results or figures. The READMEs of the folders say which input each step needs and where it comes from.

## Getting started

Scripts are run from the repository root:

```bash
Rscript NewImplementation/01-DataProcessing/ProcessData.R   # once
Rscript NewImplementation/RunComplexoFinder.R
```

See [NewImplementation/README.md](NewImplementation/README.md) for the input, the steps and the settings.
