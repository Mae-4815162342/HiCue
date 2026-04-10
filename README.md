# HiCue

**HiCue** is a command-line tool for extracting, aggregating, and visualising
chromatin interaction data from Hi-C / Micro-C experiments stored in
[Cooler](https://cooler.readthedocs.io/) (`.cool` / `.mcool`) files.

It supports pileup analysis (averaging submatrices centred on a set of genomic
positions), overlay of genomic tracks (BigWig), GFF/GTF annotation, and
multiple separation strategies (strand, region, chromosome).

---

## Table of contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Commands](#commands)
  - [extract](#extract)
  - [tracks](#tracks)
  - [compare *(experimental)*](#compare-experimental)
- [Input formats](#input-formats)
- [Output structure](#output-structure)
- [Architecture overview](#architecture-overview)
- [Known issues & fixes](#known-issues--fixes)
- [License](#license)
- [Citation](#citation)

---

## Features

- Extract Hi-C submatrices centred on single genomic loci (BED/GFF) or on
  pairs of loci (BED2D / loop anchors).
- Compute **pileup** (aggregate) matrices by median or mean across all loci in
  a group.
- Overlay **BigWig** genomic tracks on submatrix figures.
- Automatic **P(s) detrending** (distance-law normalisation) and optional
  **patch detrending** (null-model pileup subtraction).
- Separate results by **strand direction**, **chromosome**, or custom **genomic
  regions**.
- Multi-resolution support (`.mcool` files).
- Multi-threaded extraction via Python `threading` + `queue`.

---

## Requirements

- Python ≥ 3.10 (required for `match`/`case` statements)
- [cooler](https://cooler.readthedocs.io/) ≥ 0.9
- [pyBigWig](https://github.com/deeptools/pyBigWig)
- [BCBio-GFF](https://github.com/chapmanb/bcbb/tree/master/gff)
- numpy, pandas, matplotlib, scikit-learn, scipy

Install all dependencies via pip:

```bash
pip install hicue
```

---

## Installation

### From PyPI

```bash
pip install hicue
```

### From source (development)

```bash
git clone https://github.com/Mae-4815162342/HiCue.git
cd HiCue
pip install -e .
```

### Verify the installation

```bash
hicue --help
```

---

## Quick start

```bash
# Extract pileup centred on loop anchors with a 50 kb window at 1 kb resolution
hicue extract results/ anchors.bed experiment.mcool \
    --pileup \
    --windows 50000 \
    --binnings 1000

# Same with BigWig tracks overlaid
hicue tracks results/ anchors.bed experiment.mcool \
    --tracks signal.bw \
    --pileup \
    --windows 50000 \
    --binnings 1000
```

---

## Commands

### `extract`

Extract submatrices around genomic positions.

```
hicue extract OUTPUT_DIR POSITIONS COOL_FILES... [OPTIONS]
```

| Option | Type | Default | Description |
|---|---|---|---|
| `--pileup` | flag | off | Compute and display aggregate pileup figures |
| `--windows` | int… | — | Half-window size(s) in bp |
| `--binnings` | int… | — | Bin size(s) in bp (for `.mcool` multi-resolution files) |
| `--loci` | flag | off | Save individual submatrix figures |
| `--batch` | flag | off | Save batched submatrix figures (64 per page) |
| `--separate-by` | str | `""` | Separate results by `direction`, `regions`, or `chroms` |
| `--display-sense` | str | `forward` | Axis orientation (`forward` or `reverse`) |
| `--flip` | flag | off | Strand-normalise matrices (flip reverse-strand loci) |
| `--display-strand` | flag | off | Overlay transcription-direction arrows |
| `--detrending` | str | `none` | Detrending method: `none`, `ps` (distance law) or `patch` |
| `--loops` | flag | off | Treat positions as loop anchors (BED2D mode) |
| `--gff` | path | — | Annotate positions with a GFF/GTF file |
| `--threads` | int | 4 | Number of worker threads |
| `--format` | str… | `pdf` | Output figure format(s) |

### `tracks`

Same as `extract` but reads peak positions directly from a BigWig track file
(peaks are detected by thresholding the signal).

```
hicue tracks OUTPUT_DIR TRACK_FILE COOL_FILES... [OPTIONS]
```

Additional options beyond those shared with `extract`:

| Option | Type | Default | Description |
|---|---|---|---|
| `--threshold` | float | — | Absolute signal threshold for peak detection |
| `--percentage` | float | — | Keep the top N% of positions by signal |
| `--min-sep` | int | 0 | Minimum distance between retained peaks (bp) |

### `compare` *(experimental)*

Compare submatrices or pileups from two Hi-C experiments side by side and
display their log₂ ratio.  This command is currently commented out pending
refactoring.

---

## Input formats

| Format | Description |
|---|---|
| `.bed` | 3-column BED (chrom, start, end) or 6-column BED with strand |
| `.bed2d` / `.bedpe` | Paired-loci format (loop anchors): chrom1, start1, end1, chrom2, start2, end2 |
| `.gff` / `.gtf` | Gene/feature annotations for position extraction |
| `.cool` | Single-resolution Cooler Hi-C file |
| `.mcool` | Multi-resolution Cooler Hi-C file |
| `.bw` / BigWig | Continuous genomic signal tracks |

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

Intermediate files (positions, formatted pairs) can be saved with
`--save-tmp`.

---

## Architecture overview

HiCue uses a multi-threaded producer–consumer pipeline:

```
FileStreamer  ──►  Parser(s)  ──►  Annotator(s)  ──►  positions DataFrame
                                                   └──►  pairing Queue

pairing Queue  ──►  Separator(s)  ──►  PairFormater(s)  ──►  formated_pairs DataFrame

formated_pairs  ──►  Extracter  ──►  SubmatrixFormater(s)  ──►  Aggregator(s) ──►  Pileup
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
`Display`/`DisplayBatch` thread.  `asyncio.run()` creates *and destroys* an
event loop on every call.  The teardown of one loop races with `plt.close()`
called at the end of the previous coroutine, corrupting matplotlib's internal
"current figure" singleton.

**Fix (applied in `classes/AsyncDisplays.py`):**
Each worker thread now creates **one** persistent event loop
(`asyncio.new_event_loop()`) on construction and reuses it for every
rendering call via `loop.run_until_complete(coro)`.  The loop is closed in an
overridden `join()` method.

---

## License

CC BY-NC 4.0 – see [LICENSE](LICENSE) for details.

---

## Citation

If you use HiCue in your research, please cite:

```
[Citation to be added upon publication]
```
