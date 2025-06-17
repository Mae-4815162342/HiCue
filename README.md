# HiCue

A visualization tool for Hi-C datasets that enables extraction, analysis, and visualization of chromatin interaction data.

## Features

- Extract submatrices from Hi-C data around (pairs of) genomic positions
- Generate pileup visualizations of chromatin interactions
- Support for multiple file formats (BED, GFF, BigWig)
- Compare Hi-C datasets

## Installation

```bash
pip install hicue
```

For development:

```bash
git clone https://github.com/Mae-4815162342/HiCue.git
cd HiCue
pip install -e .
```

## Usage

### Extract submatrices around genomic positions

```bash
hicue extract output_dir positions.bed data.cool
```

### Generate pileup visualizations

```bash
hicue extract output_dir positions.bed data.cool --pileup
```

...

## Command Line Interface

HiCue provides several commands:

- `extract`: Extract submatrices around genomic positions
- `extract2d`: Extract submatrices for 2D genomic intervals
- `tracks`: Overlay genomic tracks on Hi-C visualizations
- `compare`: Compare two Hi-C datasets
- `annotate`: Annotate genomic positions with GFF features

Use `hicue <command> --help` for detailed options.

## License

This project is licensed under CC BY-NC 4.0. See the [LICENSE](LICENSE) file for details.

## Citation

If you use HiCue in your research, please cite:

```
...
```

