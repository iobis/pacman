# Tourmaline

Pieter Provoost, 2 March 2020

These are my notes for installing and running the Tourmaline bioinformatics pipeline. This largely follows the official documentation at https://github.com/aomlomics/tourmaline

## Installation

First clone the Tourmaline repository, then start the docker container and mount the working directory as a volume under `/data`:

```bash
git clone git@github.com:aomlomics/tourmaline.git
docker run -v "$(pwd):/data" -it --name tourmaline aomlomics/tourmaline
```

To use the container again after exit:

```bash
docker start tourmaline
docker exec -it tourmaline /bin/bash
```

I installed vim to be able to edit files from the container:

```bash
apt-get update
apt-get install vim
```

## Working with the example data

### Setting up the test data

By mounting the Tourmaline repository in `/data`, the paths in the FASTQ manifest do not need to be modified:

```bash
(qiime2-2020.8) root@4f308ff87ace:/data# head tourmaline/00-data/manifest_pe.csv
sample-id,absolute-filepath,direction
SC56.50,/data/tourmaline/00-data/fastq/Sample_134716_R1.fastq.gz,forward
SC56.50,/data/tourmaline/00-data/fastq/Sample_134716_R2.fastq.gz,reverse
SC56.22,/data/tourmaline/00-data/fastq/Sample_134717_R1.fastq.gz,forward
SC56.22,/data/tourmaline/00-data/fastq/Sample_134717_R2.fastq.gz,reverse
SC54.50,/data/tourmaline/00-data/fastq/Sample_134718_R1.fastq.gz,forward
SC54.50,/data/tourmaline/00-data/fastq/Sample_134718_R2.fastq.gz,reverse
SC54.22,/data/tourmaline/00-data/fastq/Sample_134719_R1.fastq.gz,forward
SC54.22,/data/tourmaline/00-data/fastq/Sample_134719_R2.fastq.gz,reverse
SC53.50,/data/tourmaline/00-data/fastq/Sample_134720_R1.fastq.gz,forward
```

We also need to download and link the reference database and taxonomy files:

```bash
wget https://data.qiime2.org/2020.8/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2020.8/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```

### Running Snakemake

Run Snakemake:

```bash
snakemake dada2_pe_denoise
snakemake dada2_pe_taxonomy_unfiltered
snakemake dada2_pe_diversity_unfiltered
snakemake dada2_pe_report_unfiltered

```

## Working with our own data

Todo.