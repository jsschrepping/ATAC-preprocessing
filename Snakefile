"""
Author: J. Schulte-Schrepping
Affiliation: LIMES, Uni Bonn
Aim: A simple Snakemake workflow to process paired-end stranded ATAC-Seq.
Date: 20220316
Run: snakemake --use-conda   
Docker: docker run -it --rm --name ATAC_preprocessing -v /path/to/data/:/data/ -v /path/to/index/:/data/bowtie2-indexfiles/ -v /path/to/tmp/:/tmp/ jsschrepping/bioinfo-base-image:jss_v0.0.2 /bin/bash
"""

SAMPLES=["17502", "17503", "17504", "17505", "17506", "17507", "17508", "17509", "17510", "17511", "17512", "17513", "17514", "17515", "17516", "17517", "17518", "17519", "17520", "17521", "17522", "17523", "17524", "17525", "17526", "17527", "17528", "17529", "17530", "17531", "17532", "17533", "17534", "17535", "17536", "17537", "17538", "17539", "17540"] 

READS=["R1", "R2"]

### Target rule ###

rule all:
    input:
        expand("output/qc/fastqc/Sample_{sample}_{read}_init_fastqc.html", sample=SAMPLES, read=READS),
        expand("output/trimmed/Sample_{sample}_{read}.fastq.gz", sample=SAMPLES, read=READS),
        expand("output/qc/fastqc/Sample_{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=READS),
        expand("output/mapped/Sample_{sample}.final.bam",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.final.bai",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.final.bw",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.final.flagstat",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.subsample.bam",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.subsample.bai",sample=SAMPLES),
        expand("output/mapped/Sample_{sample}.subsample.bw",sample=SAMPLES),
        expand("output/peaks/Sample_{sample}_peaks.xls",sample=SAMPLES),
        expand("output/peaks/Sample_{sample}_subsample_peaks.xls",sample=SAMPLES),
        "output/mapped/merged.bam",
        "output/peaks/merged_peaks.xls",
        "output/qc/multiqc.html"

### Initial QC ###

rule fastqc_init:
    input:
        "/data/fastq/merged_fastq/Sample_{sample}_{read}.fastq.gz"
    output:
        html="output/qc/fastqc/Sample_{sample}_{read}_init_fastqc.html",
        zip="output/qc/fastqc/Sample_{sample}_{read}_init_fastqc.zip"
    params: ""
    log:
        "output/logs/fastqc/Sample_{sample}_{read}_init.log"
    wrapper:
        "v1.3.0/bio/fastqc"

### Adapter trimming ###

rule trimmomatic_pe:
    input:
        r1="/data/fastq/merged_fastq/Sample_{sample}_R1.fastq.gz",
        r2="/data/fastq/merged_fastq/Sample_{sample}_R2.fastq.gz"
    output:
        r1="output/trimmed/Sample_{sample}_R1.fastq.gz",
        r2="output/trimmed/Sample_{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/Sample_{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="output/trimmed/Sample_{sample}_R2.unpaired.fastq.gz"
    log:
        "output/logs/trimmomatic/Sample_{sample}.log"
    params:
        trimmer=["ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36 CROP:48"],
        extra="",
        compression_level="-9"
    threads:
        16
    wrapper:
        "v1.3.0/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input:
        "output/trimmed/Sample_{sample}_{read}.fastq.gz"
    output:
        html="output/qc/fastqc/Sample_{sample}_{read}_trimmed_fastqc.html",
        zip="output/qc/fastqc/Sample_{sample}_{read}_trimmed_fastqc.zip"
    params: ""
    log:
        "output/logs/fastqc/Sample_{sample}_{read}_trimmed.log"
    wrapper:
        "v1.3.0/bio/fastqc"

### alignment ###

rule bowtie2:
    input:
        r1="output/trimmed/Sample_{sample}_R1.fastq.gz",
        r2="output/trimmed/Sample_{sample}_R2.fastq.gz"
    output:
        temp("output/mapped/Sample_{sample}.bam")
    conda:
        "envs/bowtie2.yaml"
    params:
        index="/data/bowtie2-indexfiles/genome",
        extra="-X 2000 --dovetail"
    log:
        "output/logs/bowtie2/Sample_{sample}.log"
    threads:
        8
    shell:
        "bowtie2 --threads {threads} -x {params.index} {params.extra} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -Sbh -F4 -f2 -q30 -@ {threads} - > {output}"

rule samtools_sort:
    input:
      	"output/mapped/Sample_{sample}.bam"
    output:
        temp("output/mapped/Sample_{sample}.sorted.bam")
    log:
        "output/logs/samtools/Sample_{sample}.log"
    params:
        extra="-m 4G"
    threads:
        8
    wrapper:
        "v1.3.0/bio/samtools/sort"

### ATAC processing ###

rule remove_duplicates:
    input:
        "output/mapped/Sample_{sample}.sorted.bam"
    output:
        bam="output/mapped/Sample_{sample}.dedup.bam",
        metrics="output/logs/picard_dedup/Sample_{sample}.dedup.metrics.txt"
    log:
        "output/logs/picard_dedup/Sample_{sample}.log"
    params:
        "--REMOVE_DUPLICATES true"
    wrapper:
        "v1.3.0/bio/picard/markduplicates"

rule index_dedup:
    input:
        "output/mapped/Sample_{sample}.dedup.bam"
    output:
        temp("output/mapped/Sample_{sample}.dedup.bai")
    params:
        ""
    wrapper:
        "v1.3.0/bio/samtools/index"

rule remove_offset:
    input:
        bam="output/mapped/Sample_{sample}.dedup.bam",
        bai="output/mapped/Sample_{sample}.dedup.bai"
    output:
        temp("output/mapped/Sample_{sample}.dedup.offset.bam")
    log:
        "output/logs/alignmentSieve/Sample_{sample}.log"
    conda:
        "envs/deeptools.yaml"
    threads:
        8
    shell:
        "alignmentSieve --numberOfProcessors {threads} --verbose --ATACshift --bam {input.bam} -o {output} 2> {log}"

rule sort_final:
    input:
        "output/mapped/Sample_{sample}.dedup.offset.bam"
    output:
        protected("output/mapped/Sample_{sample}.final.bam")
    log:
        "output/logs/samtools/Sample_{sample}.log"
    params:
        extra="-m 4G"
    threads:
        8
    wrapper:
        "v1.3.0/bio/samtools/sort"

rule index_final:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        protected("output/mapped/Sample_{sample}.final.bai")
    params:
        "" # optional params string
    threads:
        4
    wrapper:
        "v1.3.0/bio/samtools/index"

rule samtools_flagstat:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        "output/mapped/Sample_{sample}.final.flagstat"
    wrapper:
        "v1.3.0/bio/samtools/flagstat"

### subsample ###

rule subsample:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        protected("output/mapped/Sample_{sample}.subsample.bam")
    log:
        "output/logs/subsample/Sample_{sample}.log"
    conda:
        "envs/sambamba.yaml"
    threads:
        8
    shell:
        """
        nreads=$(samtools view -c {input})
        rate=$(echo "scale=5;2000000/$nreads" | bc)
        sambamba view -f bam -t 5 --subsampling-seed=42 -s $rate {input} | samtools sort -m 4G -@ 8 -T - > {output} 2> {log}
        """

rule index_subsample:
    input:
        "output/mapped/Sample_{sample}.subsample.bam"
    output:
        "output/mapped/Sample_{sample}.subsample.bai"
    params:
        "" # optional params string
    wrapper:
        "v1.3.0/bio/samtools/index"

### BigWig ###

rule bigwig_final:
    input:
        bam="output/mapped/Sample_{sample}.final.bam",
        bai="output/mapped/Sample_{sample}.final.bai"
    output:
        "output/mapped/Sample_{sample}.final.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"

rule bigwig_subsample:
    input:
        bam="output/mapped/Sample_{sample}.subsample.bam",
        bai="output/mapped/Sample_{sample}.subsample.bai"
    output:
        "output/mapped/Sample_{sample}.subsample.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"

### Peak Calling ###

rule macs2_full:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        xls="output/peaks/Sample_{sample}_peaks.xls"
    params:
        sample="Sample_{sample}"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/Sample_{sample}.log"
    shell:
        "macs2 callpeak -t {input} -f BAM -g hs --nomodel --extsize 50 --outdir output/peaks -n {params.sample} 2> {log}"

rule macs2_subsample:
    input:
        "output/mapped/Sample_{sample}.subsample.bam"
    output:
        xls="output/peaks/Sample_{sample}_subsample_peaks.xls"
    params:
        sample="Sample_{sample}_subsample"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/Sample_{sample}.subsample.log"
    shell:
        "macs2 callpeak -t {input} -f BAM -g hs --nomodel --extsize 50 --outdir output/peaks -n {params.sample} 2> {log}"


### Merge bam files

rule samtools_merge:
    input:
        expand("output/mapped/Sample_{sample}.final.bam", sample=SAMPLES),
    output:
        "output/mapped/merged.bam",
    log:
        "output/logs/samtools/merge.log",
    params:
        extra="",  # optional additional parameters as string
    threads: 8
    wrapper:
        "v1.3.1/bio/samtools/merge"

### Call peaks on merged bam file

rule macs2_merged:
    input:
        "output/mapped/merged.bam"
    output:
        xls="output/peaks/merged_peaks.xls"
    params:
        sample="merged"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/merged.log"
    shell:
        "macs2 callpeak -t {input} -f BAM -g hs --outdir output/peaks -n {params.sample} 2> {log}"

### Complete QC ###

rule multiqc:
    input:
        expand("output/qc/fastqc/Sample_{sample}_{read}_init_fastqc.zip",zip,sample=SAMPLES,read=READS),
        expand("output/qc/fastqc/Sample_{sample}_{read}_trimmed_fastqc.zip",zip,sample=SAMPLES,read=READS),
        expand("output/logs/bowtie2/Sample_{sample}.log", sample=SAMPLES),
        expand("output/peaks/Sample_{sample}_peaks.xls", sample=SAMPLES),
        expand("output/peaks/Sample_{sample}_subsample_peaks.xls", sample=SAMPLES)
    output:
        "output/qc/multiqc.html"
    log:
        "output/logs/multiqc.log"
    params:
        ""
    wrapper:
        "v1.3.0/bio/multiqc"
