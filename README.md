# Demultiplexing of NanoStrand-seq

NanoStrand-seq utilized the dual-barcode combination to split single-cell reads from mixed reads and determine the read direction.

## Installation

    # Added following line to your '.bashrc' file:
    export PATH="`pwd`/scripts:${PATH}"

## Preparation

Firstly, we should prepare the following 4 files in advance:

No.|File|Description
:-|:-|:-
1|reads.fastq.gz|Mixed sequencing reads in FASTQ format.
2|barcodes1.fasta|Sequences of the 1st barcode in FASTA format.
3|barcodes2.fasta|Sequences of the 2nd barcode in FASTA format.
4|barcode_config.tsv|A tab-delimited file that describe the cell and corresponding barcode combination.

## Running

We provied an Snakemake workflow template for demultiplexing single-cell reads in `snakemake.smk`.

The required files were deposited in the `./data` directory.

You can run the workflow by following command:

    snakemake -s snakemake.smk -j 4

## How to determine threshold?



