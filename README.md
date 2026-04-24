# HiCue

**HiCue** is a command-line tool for extracting, aggregating, and visualising
chromatin interaction data from Hi-C / Micro-C experiments stored in
[Cooler](https://cooler.readthedocs.io/) (`.cool` / `.mcool`) files.

It supports pileup analysis (averaging submatrices centred on a set of genomic
positions or regions), overlay of genomic tracks (BigWig), GFF/GTF annotation,
and multiple separation strategies (strand, region, chromosome).

---

## Table of contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Commands](#commands)
  - [extract](#extract)
  - [tracks](#tracks)
  - [regions](#regions)
- [Input formats](#input-formats)
- [Passing Cooler files](#passing-cooler-files)
- [Output structure](#output-structure)
- [Examples](#examples)
- [Architecture overview](#architecture-overview)
- [Known issues & fixes](#known-issues--fixes)
- [License](#license)
- [Citation](#citation)

---

## Features

- Extract Hi-C submatrices centred on single genomic loci (BED/GFF) or on
  pairs of loci (BED2D / loop anchors).
- Compute **pileup** (aggregate) matrices by median, mean, or sum across all
  loci in a group.
- Extract and aggregate submatrices of **variable genomic size** using the
  `regions` command, which rescales each submatrix to a common pixel
  dimension before aggregation.
- Overlay **BigWig** genomic tracks on submatrix figures.
- Automatic **P(s) detrending** (distance-law normalisation) and optional
  **patch detrending** (null-model pileup subtraction).
- Separate results by **strand direction**, **chromosome**, or custom
  **genomic regions**.
- Multi-resolution support (`.mcool` files).
- Pass multiple Cooler files at once — directly on the command line, as a
  comma-separated list, or via a **plain-text file** listing one path per line.
- Multi-threaded extraction via Python `threading` + `queue`.

---

## Requirements

- Python ≥ 3.10
- [cooler](https://cooler.readthedocs.io/) ≥ 0.10
- [pyBigWig](https://github.com/deeptools/pyBigWig)
- [BCBio-GFF](https://github.com/chapmanb/bcbb/tree/master/gff)
- numpy, pandas, matplotlib, chromosight, click

All dependencies are installed automatically when you install HiCue via pip.

---

## Installation

### Option A — pip (recommended)

```bash
pip install hicue
```

It is strongly recommended to install inside a dedicated environment:

```bash
python -m venv hicue-env
source hicue-env/bin/activate   # Linux / macOS
# hicue-env\Scripts\activate   # Windows
pip install hicue
```

### Option B — conda environment

```bash
conda create -n hicue python=3.11
conda activate hicue
pip install hicue
```

> **Why pip inside conda?**  Some HiCue dependencies (e.g. `pyBigWig`) are not
> available on the default conda channels for all platforms. Using pip inside a
> conda environment gives you the isolation of conda with the full package
> availability of PyPI.

### Option C — from source (development)

```bash
git clone https://github.com/Mae-4815162342/HiCue.git
cd HiCue
pip install -e .[dev]
```

### Verify

```bash
hicue --help
```

---

## Quick start

```bash
# 1 — Install
pip install hicue

# 2 — Pileup centred on loop anchors, 50 kb window, 1 kb resolution
hicue extract results/ anchors.bed experiment.mcool \
    --windows 50000 \
    --binnings 1000

# 3 — Same analysis on several experiments at once (text-file input)
hicue extract results/ anchors.bed experiments.txt \
    --windows 50000 \
    --binnings 1000

# 4 — Variable-size region pileup (e.g. TADs)
hicue regions results/ tads.bed experiment.mcool \
    --expected_sizes 51 \
    --padding 0.5 \
    --binnings 5000
```

---

## Commands

### `extract`

Extract submatrices around fixed-size genomic positions and compute a pileup.

```
hicue extract OUTPUT_DIR POSITIONS COOL_FILES [OPTIONS]
```

| Option | Type | Default | Description |
|---|---|---|---|
| `--pileup / --no-pileup` | flag | on | Compute and display aggregate pileup figures |
| `--loci / --no-loci` | flag | off | Save individual submatrix figures |
| `--batch / --no-batch` | flag | off | Save batched submatrix figures (64 per page) |
| `-w / --windows` | int,… | 30000 | Half-window size(s) in bp |
| `-b / --binnings` | int,… | 1000 | Bin size(s) in bp (`.mcool` only) |
| `-d / --detrending` | choice | `none` | `none`, `ps` (distance law), or `patch` (null-model subtraction) |
| `-m / --method` | choice | `median` | Aggregation method: `median`, `mean`, or `sum` |
| `-f / --flip` | flag | off | Strand-normalise matrices (flip reverse-strand loci) |
| `-r / --raw` | flag | off | Use raw (unbalanced) contact counts |
| `--loops` | flag | off | Treat positions as loop anchors (BED2D mode) |
| `--trans` | flag | off | Include trans-chromosomal contacts in loop mode |
| `--min_dist` | int | 30000 | Minimum distance in bp between paired loci |
| `--diag_mask` | int | 0 | Mask the diagonal up to this distance (bp) |
| `--separate_by` | str,… | — | Separate results by `direction`, `regions`, or `chroms` |
| `--separation_regions` | path | — | CSV file defining regions for `--separate_by regions` |
| `--gff` | path | — | Annotate positions with a GFF/GTF file |
| `--tracks` | path | — | Annotate positions with a BigWig file |
| `--center` | choice | `start` | Window anchor: `start`, `center`, or `end` of each feature |
| `--display_sense` | choice | `forward` | Axis orientation: `forward` or `reverse` |
| `--display_strand` | flag | off | Overlay transcription-direction arrows on figures |
| `--cmap_color` | str | `seismic` | Matplotlib colormap for pileup figures |
| `--cmap_limits` | float float | — | Fixed min/max for the pileup colormap |
| `--indiv_cmap_color` | str | `afmhot_r` | Colormap for individual submatrix figures |
| `--indiv_cmap_limits` | float float | — | Fixed min/max for individual figures |
| `--format` | str,… | `pdf` | Output figure format(s) (e.g. `pdf,png`) |
| `--threads` | int | 8 | Number of worker threads |
| `--nb_pos` | int | 2 | Random positions per locus for patch detrending |
| `--rand_max_dist` | int | 100000 | Max distance (bp) for random position selection |
| `--random_jitter` | int | 0 | Allowed jitter (bp) on random pair distances |
| `--random_path` | str | — | Path prefix to pre-computed random positions |
| `--circulars` | str,… | `none` | Chromosomes to treat as circular |
| `--ps_all_chrom` | bool | True | Use all cis contacts for P(s) estimation |
| `--contact_range` | int int int | — | MIN MAX STEP (bp) for distance-based separation |
| `--overlap` | choice | `strict` | Position-interval overlap: `strict` or `flex` |
| `--record_type` | choice | — | GFF record type to select (e.g. `gene`) |
| `--track_unit` | str | `""` | Label for the BigWig signal axis |
| `--save_tmp` | flag | off | Save intermediate CSV files to OUTPUT_DIR |

### `tracks`

Equivalent to `extract` but derives genomic positions directly from peaks
detected in a BigWig signal file.

```
hicue tracks OUTPUT_DIR TRACK_FILE COOL_FILES [OPTIONS]
```

All options from `extract` are supported. Additional options:

| Option | Type | Default | Description |
|---|---|---|---|
| `-t / --threshold` | choice+float | — | `min VALUE` keeps positions above VALUE; `max VALUE` keeps positions below |
| `-p / --percentage` | choice+int | — | `high N` keeps the top N%, `low N` keeps the bottom N% of positions by signal |
| `--min_sep` | int | 1000 | Minimum distance in bp between two retained peaks |
| `--positions` | path | — | Restrict peak selection to positions in this BED/GFF file |
| `--gff_type` | str | `""` | Feature type to select when `--positions` is a GFF file |

### `regions`

Extract and aggregate submatrices of **variable genomic size**, resizing each
one to a common pixel dimension before computing the pileup. Designed for
features with different lengths (TADs, genes, compartments).

```
hicue regions OUTPUT_DIR POSITIONS COOL_FILES [OPTIONS]
```

All options from `extract` are supported. Additional options:

| Option | Type | Default | Description |
|---|---|---|---|
| `-p / --padding` | float | 1.0 | Padding ratio added on each side of the region (e.g. 0.5 adds half the region size) |
| `-s / --min_region_size` | int | 20000 | Minimum region size in bp; smaller regions are skipped |
| `-e / --expected_sizes` | int,… | 51 | Target pixel dimension(s) for resizing (e.g. `20,51`) |

> **Note:** Small regions combined with a large bin size can introduce strong
> biases in the pileup. Use `--min_region_size` to filter them out and prefer
> a bin size that yields at least ~5 bins per region.

---

## Input formats

| Format | Extension(s) | Notes |
|---|---|---|
| BED | `.bed` | 3-column minimum. Column 4 = name, column 6 = strand (`+`/`-`). Strand required for `--flip` and `--separate_by direction`. |
| BED2D / BEDPE | `.bed2d`, `.bedpe` | 6-column: chrom1, start1, end1, chrom2, start2, end2. Used with `--loops`. |
| GFF / GTF | `.gff`, `.gtf` | Gene/feature annotations for position extraction or labelling. |
| Cooler (single) | `.cool` | Single-resolution matrix. `--binnings` is ignored. |
| Cooler (multi) | `.mcool` | Multi-resolution matrix. Select resolution(s) with `--binnings`. |
| BigWig | `.bw` | Continuous signal track. Primary input for `tracks`, or annotation overlay in `extract`/`regions`. |
| Regions CSV | `.csv` | Comma-separated: `Id,Chromosome,Start,End`. Used with `--separate_by regions`. Multiple rows with the same `Id` define a discontinuous interval. |

---

## Passing Cooler files

The `COOL_FILES` argument is flexible and accepts three forms:

**1. A single file directly:**
```bash
hicue extract results/ anchors.bed experiment.mcool --binnings 1000
```

**2. A comma-separated list of files:**
```bash
hicue extract results/ anchors.bed control.mcool,treated.mcool --binnings 1000
```

**3. A plain-text file listing one Cooler path per line:**
```
# experiments.txt
/data/project/control.mcool
/data/project/treated.mcool
/data/project/recovery.mcool
```
```bash
hicue extract results/ anchors.bed experiments.txt --binnings 1000
```

All three forms work identically with `extract`, `tracks`, and `regions`.
Each Cooler file produces its own sub-directory inside `OUTPUT_DIR`.

---

## Output structure

```
OUTPUT_DIR/
└── {cool_name}/
    └── {sep_id}/
        └── binning_{binning}/
            ├── individual_{window}kb_window/
            │   ├── GeneName.pdf          ← per-locus submatrix
            │   └── …
            ├── batched_{window}kb_window/
            │   ├── batch#1.pdf           ← 64-panel batch figure
            │   ├── batch#1_references.csv
            │   └── …
            └── pileup_{window}kb_window.pdf
```

A timestamped log file (`YYYYMMDD_HHMMSS_log.txt`) recording the exact command
and all option values is written to `OUTPUT_DIR` for every run.

Intermediate CSV files (positions, formatted pairs, random positions for patch
detrending) can be saved with `--save_tmp` and reused in subsequent runs via
`--random_path`.

---

## Examples

### Fixed-window pileup on gene TSS

```bash
hicue extract results/tss/ genes.bed experiment.mcool \
    --windows 50000 \
    --binnings 1000 \
    --center start
```

### Strand-aware pileup with individual figures

```bash
hicue extract results/tss_stranded/ genes.bed experiment.mcool \
    --windows 50000 \
    --binnings 1000 \
    --center start \
    --flip \
    --loci
```

### Loop pileup with P(s) detrending, multiple windows

```bash
hicue extract results/loops/ loops.bed2d control.mcool treated.mcool \
    --loops \
    --detrending ps \
    --windows 50000,100000 \
    --binnings 1000 \
    --format pdf,png
```

### Patch detrending — generate then reuse random positions

```bash
# First run: compute and save random control positions
hicue extract results/ctrl/ loops.bed2d control.mcool \
    --loops \
    --detrending patch \
    --nb_pos 3 \
    --save_tmp

# Second run: reuse the same random positions for a paired comparison
hicue extract results/treated/ loops.bed2d treated.mcool \
    --loops \
    --detrending patch \
    --nb_pos 3 \
    --random_path results/ctrl/loops
```

### Pileup from BigWig peaks

```bash
# Keep the top 20 % of ChIP-seq signal bins, minimum 5 kb between peaks
hicue tracks results/chip/ H3K27ac.bw experiment.mcool \
    --percentage high 20 \
    --min_sep 5000 \
    --windows 50000 \
    --binnings 1000
```

### Variable-size region pileup (TADs)

```bash
hicue regions results/tads/ tads.bed experiment.mcool \
    --expected_sizes 51 \
    --padding 0.5 \
    --binnings 5000
```

### Multiple resolutions from a file list

```bash
hicue regions results/tads_multi/ tads.bed experiments.txt \
    --expected_sizes 20,51 \
    --binnings 5000,10000 \
    --padding 1.0
```

### Separating results by chromosome

```bash
hicue extract results/by_chrom/ genes.bed experiment.mcool \
    --windows 50000 \
    --binnings 5000 \
    --separate_by chroms
```

---

## Architecture overview

HiCue uses a multi-threaded producer–consumer pipeline:

```
FileStreamer  ──►  Parser(s)  ──►  Annotator(s)  ──►  positions DataFrame
                                                   └──►  pairing Queue

pairing Queue  ──►  Separator(s)  ──►  PairFormater(s)  ──►  formated_pairs DataFrame

formated_pairs  ──►  Extracter  ──►  SubmatrixFormater(s)  ──►  Aggregator(s)  ──►  Pileup
                                  └──►  DisplayBatch / Display  (async rendering)

Pileup Queue  ──►  Display (pileup rendering)
```

All inter-stage communication uses `queue.Queue` with a `"DONE"` sentinel.
Display workers (in `classes/AsyncDisplays.py`) run `async` matplotlib
functions on a dedicated per-thread event loop.

---

## Known issues & fixes

### Blank / corrupted figures with async display workers

**Symptom:** Output PDF/PNG files are blank, show the wrong data, or the
process crashes with a `RuntimeError` from the Agg renderer.

**Root cause:** The original code used `asyncio.run()` inside each
`Display`/`DisplayBatch` thread. `asyncio.run()` creates *and destroys* an
event loop on every call. The teardown of one loop races with `plt.close()`
called at the end of the previous coroutine, corrupting matplotlib's internal
"current figure" singleton.

**Fix (applied in `classes/AsyncDisplays.py`):**
Each worker thread now creates **one** persistent event loop
(`asyncio.new_event_loop()`) on construction and reuses it for every
rendering call via `loop.run_until_complete(coro)`. The loop is closed in an
overridden `join()` method.

---

## License

CC BY-NC 4.0 – see [LICENSE](LICENSE) for details.

---
## Acknowledgments

A. Cournac for project supervision and tests.
J. Serizay for primary documentation and test implementation.
M. Perrot for provided data.

## Citation

If you use HiCue in your research, please cite:

```
[Citation to be added upon publication]
```