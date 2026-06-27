```
█▀▀█ █░░█ █▀▀█ █▀▀▀ █▀▀ 　 █▀▀ █▀▀█ █▀▀▄ ▀▀█▀▀ ░▀░ █▀▀▀ 　 █▀▀█ █▀▀▄ █▀▀▄ █▀▀█ ▀▀█▀▀ █▀▀█ ▀▀█▀▀ █▀▀█ █▀▀█ 
█░░█ █▀▀█ █▄▄█ █░▀█ █▀▀ 　 █░░ █░░█ █░░█ ░░█░░ ▀█▀ █░▀█ 　 █▄▄█ █░░█ █░░█ █░░█ ░░█░░ █▄▄█ ░░█░░ █░░█ █▄▄▀ 
█▀▀▀ ▀░░▀ ▀░░▀ ▀▀▀▀ ▀▀▀ 　 ▀▀▀ ▀▀▀▀ ▀░░▀ ░░▀░░ ▀▀▀ ▀▀▀▀ 　 ▀░░▀ ▀░░▀ ▀░░▀ ▀▀▀▀ ░░▀░░ ▀░░▀ ░░▀░░ ▀▀▀▀ ▀░▀▀
```

Annotates genes on putative phage contigs with a database of hidden Markov models (HMMs) based on [PHROGs](https://phrogs.lmge.uca.fr/). This tool was built to support visual confirmation of predictions made by [Jaeger](https://github.com/Yasas1994/Jaeger).

The pipeline is managed by [Snakemake](https://snakemake.github.io/), uses [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) for gene calling, [pyhmmer](https://github.com/althonos/pyhmmer) for HMM search, and [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) for tRNA prediction.

## Requirements

- Python >= 3.11
- Conda (recommended for installing tRNAscan-SE)

## Installation

1) Clone the repository:
   ```bash
   git clone https://github.com/Yasas1994/phage_contig_annotator.git
   cd phage_contig_annotator
   ```

2) Create and activate the Conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate phage_contig_annot
   ```

3) Install the package in editable mode:
   ```bash
   pip install -e ".[dev]"
   ```

4) Download the annotation database:
   ```bash
   phage_contig_annotator download-db
   ```

## Usage

Activate the environment and run the full pipeline:

```bash
conda activate phage_contig_annot
phage_contig_annotator runall --input test/bin.460.fna --output output_dir --cpus 10
```

To skip tRNA prediction:

```bash
phage_contig_annotator runall --input test/bin.460.fna --output output_dir --cpus 10 --skip-trna
```

To run on an existing protein FASTA:

```bash
phage_contig_annotator runall --input proteins.faa --output output_dir --type proteins --cpus 10
```

To preview the Snakemake execution plan without running:

```bash
phage_contig_annotator runall --input test/bin.460.fna --output output_dir --cpus 10 --dry-run
```

## Development

Run the test suite:

```bash
pytest tests/ -v
```

## Examples

![image](https://user-images.githubusercontent.com/34155351/201523981-1f5f2fd0-5177-417e-aea0-752430c7eb80.png)
![image](https://user-images.githubusercontent.com/34155351/201524028-13f35b46-83fb-4c50-af81-3701ef992e48.png)
![image](https://user-images.githubusercontent.com/34155351/201524110-c4919fd4-aa34-449b-9ba1-7081a278d6af.png)
