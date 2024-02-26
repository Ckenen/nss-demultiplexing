#!/usr/bin/env snakemake

runs = ["20220708_GM12878"]
cells = [line.strip().split()[0] for line in open("data/20220708_GM12878.barcode_config.tsv")]
run_cells = ["20220708_GM12878/%s" % c for c in cells]
# print(run_cells)

rule all:
    input:
        expand("results/fbilr/{run}.matrix.gz", run=runs),
        expand("results/fbilr/{run}.stats.tsv.gz", run=runs),
        expand("results/splitted/{run}", run=runs),
        expand("results/combined/{run_cell}.fastq.gz", run_cell=run_cells),
        expand("results/trimmed/{run_cell}", run_cell=run_cells),

rule fbilr:
    input:
        fq = "data/{run}.10k.fastq.gz",
        fa1 = "data/{run}.1st.fasta",
        fa2 = "data/{run}.2nd.fasta"
    output:
        tsv = "results/fbilr/{run}.matrix.gz"
    log:
        "results/fbilr/{run}.log"
    threads:
        4
    shell:
        """
        fbilr -t {threads} -w 200 -b {input.fa1},{input.fa2} {input.fq} 2> {log} \
            | gzip -c > {output.tsv}
        """

rule stats:
    input:
        mtx = "results/fbilr/{run}.matrix.gz"
    output:
        tsv = "results/fbilr/{run}.stats.tsv.gz"
    shell:
        """
        zcat {input.mtx} | awk '$2>=400' | awk -v OFS=',' '{{print $3,$4,$5,$8,$9,$10,$11,$14}}' \
            | sort | uniq -c | awk -v OFS=',' '{{print $2,$1}}' \
            | sed 's/,/\\t/g' | gzip -c > {output.tsv}
        """

rule split:
    input:
        fq = "data/{run}.10k.fastq.gz",
        tsv = "results/fbilr/{run}.matrix.gz",
        cfg = "data/{run}.barcode_config.tsv"
    output:
        out = directory("results/splitted/{run}")
    log:
        "results/splitted/{run}.log"
    shell:
        """
        ./scripts/nss_split_reads.py -e 6 -l 400 {input.fq} {input.tsv} {input.cfg} {output.out} &> {log}
        """

rule combine:
    input:
        fqs = "results/splitted/{run}"
    output:
        fq = "results/combined/{run}/{cell}.fastq.gz"
    threads:
        4
    shell:
        """
        fq1={input.fqs}/fastqs/{wildcards.cell}_F.fastq
        fq2={input.fqs}/fastqs/{wildcards.cell}_R.fastq
        ( cat $fq1; cat $fq2 | ./scripts/reverse_fastq.py ) | pigz -p {threads} -c > {output.fq}
        """

rule trim:
    input:
        fq = "results/combined/{run}/{cell}.fastq.gz"
    output:
        out = directory("results/trimmed/{run}/{cell}")
    log:
        "results/trimmed/{run}/{cell}.log"
    shell:
        """
        ./scripts/nss_trim_reads.py {input.fq} {output.out} &> {log}
        """
