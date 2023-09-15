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

Here, we show an example of the above file in the `data` directory.

The `reads.fastq.gz` is located at `data/20220708_GM12878.10k.fastq.gz`, and the content is:

    @acbfb367-2b0d-48de-913b-72e4730574c1 runid=fa2a5215cdbe4ea7685499d65557d30979bd24cd read=9 ch=90 start_time=2022-07-12T14:34:22.236730+08:00 flow_cell_id=PAK14180 protocol_group_id=FKDN220228112-1A sample_id=FOND220228112-1A parent_read_id=acbfb367-2b0d-48de-913b-72e4730574c1
    GTACTTCGTTCAGTTACGTATTGCTACACTCTTTCCCTACACGACGCTCTTCAGGGAACAAACCAAGTTACGTTCGTCGGCAGCGTCAGATGTGTATAAGGAGACAGGAGGCTGAGGCATGAAAATCACTTGAACCAGGGAGGTGGGGGTTACAGTGAGCTGAAATCATGCCACTGCACTCCAGCCTGGGTGACAGAGAGAGACTGTCTCAAAAAAAAAAAAA
    +
    &'((23322:=<;>==<9:90.))''84346?ADD@>==@77<=;<<;:<<53220589<@==<<657>=<-199876,++-973335::::><=<<C65544=:7867;:89888@??@B@++++97888@::????>ABC=876****65658:<=>B=100;;<??=<==ACCAC@=@>?=<>>A@?9899?>?H{{{{{A=8876323A>DDKC@<73.

    ...

The `barcodes1.fa` is located at `data/20220708_GM12878.p7.fasta`, and the content is:

    >Bar7
    GTGTTACCGTGGGAATGAATCCTT
    >Bar8
    TTCAGGGAACAAACCAAGTTACGT
    >Bar19
    GTTCCTCGTGCAGTGTCAAGAGAT
    >Bar20
    TTGCGTCCTGTTACGAGAACTCAT
    >Bar21
    GAGCCTCTCATTGTCCGTTCTCTA
    >Bar22
    ACCACTGCCATGTATCAAAGTACG
    >Bar23
    CTTACTACCCAGTGAACCTCCTCG
    >Bar26
    CATACAGCGACTACGCATTCTCAT
    >Bar27
    CGACGGTTAGATTCACCTCTTACA
    >Bar28
    TGAAACCTAAGAAGGCACCGTATC
    >Bar29
    CTAGACACCTTGGGTTGACAGACC
    >Bar30
    TCAGTGAGGATCTACTTCGACCCA

The `barcodes2.fa` is located at `data/20220708_GM12878.p5.fasta`, and the content is:

    >Bar9
    AACTAGGCACAGCGAGTCTTGGTT
    >Bar11
    GTTTCATCTATCGGAGGGAATGGA
    >Bar12
    CAGGTAGAAAGAAGCAGAATCGGA
    >Bar13
    AGAACGACTTCCATACTCGTGTGA
    >Bar14
    AACGAGTCTCTTGGGACCCATAGA
    >Bar15
    AGGTCTACCTCGCTAACACCACTG
    >Bar16
    CGTCAACTGACAGTGGTTCGTACT
    >Bar18
    CCAAACCCAACAACCTAGATAGGC

The `barcode_config.tsv` is located at `data/20220708_GM12878.barcode_config.tsv`, and the content is:

    20220708_GM12878.sc001	Bar19	Bar9
    20220708_GM12878.sc002	Bar19	Bar11
    20220708_GM12878.sc003	Bar19	Bar12
    20220708_GM12878.sc004	Bar19	Bar13
    20220708_GM12878.sc005	Bar19	Bar14
    20220708_GM12878.sc006	Bar19	Bar15
    20220708_GM12878.sc007	Bar19	Bar16
    20220708_GM12878.sc008	Bar19	Bar18
    20220708_GM12878.sc009	Bar20	Bar9
    20220708_GM12878.sc010	Bar20	Bar11

    ...

## Running

Firstly, we find the barcodes in the reads by using FBILR:

    fbilr -t 4 -w 200 -b data/20220708_GM12878.p7.fasta,data/20220708_GM12878.p5.fasta data/20220708_GM12878.10k.fastq.gz | gzip -c > data/20220708_GM12878.matrix.tsv.gz

The content of `data/20220708_GM12878.matrix.tsv.gz` is:

    ...

    94a0332a-be37-4387-83c4-81e32165fb10	3291	Bar16	R	T	3227	3252	2	Bar7	F	T	3259	3278
    1b36b5f5-2d05-4a89-9499-c800d65509d3	2277	Bar16	F	H	52	79	6	Bar7	R	T	2212	2230
    e8f9414a-23f8-463c-a56f-cf8651c8fa9d	1572	Bar18	R	T	1521	1543	5	Bar27	F	H	44	68 0
    86b3704d-e0af-4a97-a649-3b82a9e69e52	1387	Bar13	F	H	47	71	0	Bar27	R	T	1329	1350
    ee573f16-03d6-40c8-b977-3703b956a481	3258	Bar12	F	H	55	80	1	Bar20	R	T	3201	3220
    d1039ac0-33d9-4ae0-86d1-60d3caf3768d	3611	Bar18	R	T	3551	3574	2	Bar8	F	H	48	72 0
    ec50b3a4-8f59-4308-b006-1a3aed9b2ee8	2736	Bar18	R	T	2674	2701	3	Bar23	F	T	2617	2647
    822682d5-8234-41dc-a88e-a239931c8e8a	2700	Bar14	R	H	28	51	7	Bar23	R	T	2637	2661
    b101d65d-08ff-4de0-ab1f-d0b49a83054a	2201	Bar12	F	H	47	71	3	Bar21	R	T	2139	2161
    81ad8ab0-5517-48ea-9ad6-ac60d9910518	2501	Bar13	F	H	52	76	0	Bar28	R	T	2442	2461

    ...

Then, we split the single-cell reads by the following command:

    ./scripts/nss_split_reads.py -f data/20220708_GM12878.10k.fastq.gz -m data/20220708_GM12878.matrix.tsv.gz -b data/20220708_GM12878.barcode_config.tsv -o data/splitted.d -r data/splitted.reads.tsv -e 5 -l 400 > data/splitted.log

    # For one cell
    cat data/splitted.d/20220708_GM12878.sc001_F.fastq > data/combined/20220708_GM12878.sc001.fastq
    cat data/splitted.d/20220708_GM12878.sc001_R.fastq | ./scripts/reverse_fastq.py >> data/combined/20220708_GM12878.sc001.fastq
    pigz -p 8 data/combined/20220708_GM12878.sc001.fastq

    # For all cells
    for fq in data/splitted.d/*_F.fastq; do
        cell=`basename $fq _F.fastq`
        cat data/splitted.d/${cell}_F.fastq > data/combined/${cell}.fastq
        cat data/splitted.d/${cell}_R.fastq | ./scripts/reverse_fastq.py >> data/combined/${cell}.fastq
        pigz -p 8 data/combined/${cell}.fastq
    done


Finally, we remove the linker sequences and potential chimeric reads by the following command: 

    ./scripts/nss_trim_reads.py data/combined/20220708_GM12878.sc001.fastq.gz data/trimmed/20220708_GM12878.sc001.fastq.gz > data/trimmed/20220708_GM12878.sc001.log

## How to determine threshold?



