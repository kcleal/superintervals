# Reproducible benchmark

Run from the repository root with Python 3.10+, a C++17 compiler and Rust
(CI pins Python 3.12, GCC 13 and Rust 1.87.0):

```sh
python3 test/benchmarking/run_benchmarks.py --output benchmark-results
```

The default runs all eight datasets in sorted and seeded shuffled order. Real
inputs retain a deterministic 25% sample across each file, after filtering to
nonempty chr1 intervals. Reciprocal comparisons use identical samples.
Mitochondrial reads are restricted to starts in the first 1,000 bp; interval ends
are preserved. Use `--sample-fraction 0.5` or `1.0` for larger inputs and
`--mito-window 0` for the whole mitochondrial chromosome. Sampling both sides
reduces coverage and output volume; parameters and exact inputs are recorded.
`--smoke` caps each sampled input at 1,000 records using reservoir sampling;
`--max-records` offers another cap. Smoke timings are validation only.

The fixed-width simulation defaults to 250,000 references and 25,000 queries,
all 150 bp, with independent uniform starts in a 3.75 Mb window (10× reference
coverage). Counts and window scale with the sampling fraction. Seeds are
`--seed` and `--seed + 1`.

The runner downloads and verifies the release archive, builds SuperIntervals
C++/Rust, IITree, IntervalTree, NCLS, COITrees, rust-lapper and sif-itree, then checks results
against independent hit-count and payload-checksum oracles. BED coordinates are
half-open; inclusive APIs receive end−1. Find fills a reusable vector; stream
consumes borrowed values without making a list; count uses specialised APIs. sif-itree uses its serial callback API for Find
and Stream; Count is unavailable because it has no specialised counting API.
Compare tools within each operation. Reciprocal directions appear together.

One output folder contains `report.html`, Markdown, embedded timing heatmaps,
raw and summary CSVs, `runtimes.csv`, metadata, logs, source snapshots and exact
prepared inputs. It must be new. Query timings exclude I/O; construction is
reported separately. Medians and Q1–Q3 summarise repeated measurements.
Elapsed process times include parsing, construction, all supported operations
and validation; they are not comparable query timings. Console progress counts
completed tool processes across all selected datasets, both orders, warmups and
measured repeats. The percentage measures completed jobs, not elapsed time. Measure one repeat with
`--repeats 1 --warmups 0`; regenerate reports using `--report-only`.

## GitHub Actions and Pages

The workflow runs on relevant pushes to main and can be triggered manually.
It defaults to all datasets at 25%, five measured repetitions and one warmup,
on both Ubuntu x86-64 and Arm64. Manual inputs select architecture, fraction,
repetitions and full/compact/smoke suites. Compact selects rna–annotation,
genes–DB53, ONT–DB53 and fixed-width. CPU model, architecture, native compiler
flags, toolchains and input hashes are recorded; hosted CPU models can vary.

Enable **Settings → Pages → Source: GitHub Actions** once. Successful non-smoke
runs on the default branch publish one page with architecture tabs. The Actions
summary links the page. CSVs and metadata are accessible there without an
artifact download; timing history retains the latest 30 runs in `history.json`.
Source snapshots and prepared inputs remain downloadable Actions artifacts for
30 days. Results are generated, never committed, and require no separate branch
or repository. Failed benchmark jobs leave the last successful page intact.

Dataset: release v0.2.0 `data.tar.gz`, SHA256
`c98159506689c412a909505e79576035b42ab50a4930d6c693ef9517756039de`.
C/C++ and NCLS use tracked sources from `test/3rd-party`; their hashes are
recorded. Cargo downloads COITrees 0.4.0, rust-lapper 1.3.0 and sif-itree 0.4.1
from crates.io.
Exact versions and transitive dependency checksums are pinned in the driver
`Cargo.lock`; builds use `--locked`. No local Rust third-party checkout or
submodule is required. Resolved package versions and sources are recorded in
metadata and the HTML tool table. Unversioned C/C++ sources show their last
Git modification date (not an upstream release date). CI fetches full history
for these dates. Cargo's full dependency graph is recorded in `cargo-metadata.json`.
To upgrade a Rust dependency, change its pinned version in the generated driver
manifest and update `test/benchmarking/Cargo.lock` together.
