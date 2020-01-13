"""
Author: J. Schulte-Schrepping
Affiliation: LIMES, Uni Bonn
Aim: A simple Snakemake workflow to process paired-end stranded ATAC-Seq.
Date: 20190528
Run: snakemake --use-conda   
Docker: docker run -it --rm --name ATAC_preprocessing -v /path/to/data/:/data/ -v /path/to/index/:/data/bowtie2-indexfiles/ -v /path/to/tmp/:/tmp/ jsschrepping/bioinfo-base-image:jss_v0.0.2 /bin/bash
"""

SAMPLES=["1111", "2222", "3333"] #exchange sample IDs according to your data
READS=["R1", "R2"]

### Target rule ###

rule all:
    input:
        expand("output/qc/fastqc/Sample_{sample}_{read}_fastqc.html", sample=SAMPLES, read=READS),
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
        "output/qc/multiqc_raw.html",
        "output/qc/multiqc.html"

### raw fastq ###

rule fastqc_raw:
    input:
        "/data/fastq/merged_fastq/Sample_{sample}_{read}.fastq.gz"
    output:
        html="output/qc/fastqc/Sample_{sample}_{read}_fastqc.html",
        zip="output/qc/fastqc/Sample_{sample}_{read}_fastqc.zip"
    params: ""
    log:
        "output/logs/fastqc/Sample_{sample}_{read}.log"
    wrapper:
        "0.34.0/bio/fastqc"

rule multiqc_raw:
    input:
        expand("output/qc/fastqc/Sample_{sample}_{read}_fastqc.zip",zip,sample=SAMPLES,read=READS)
    output:
        "output/qc/multiqc_raw.html"
    log:
        "output/logs/multiqc_raw.log"
    params:
        ""
    wrapper:
        "0.34.0/bio/multiqc"

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
        trimmer=["ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36"],
        extra="",
        compression_level="-9"
    threads:
        32
    wrapper:
        "0.35.0/bio/trimmomatic/pe"

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
        "0.34.0/bio/fastqc"

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
        "-m 4G"
    threads:
        8
    wrapper:
        "0.35.0/bio/samtools/sort"

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
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.35.0/bio/picard/markduplicates"

rule index_dedup:
    input:
        "output/mapped/Sample_{sample}.dedup.bam"
    output:
        temp("output/mapped/Sample_{sample}.dedup.bai")
    params:
        ""
    wrapper:
        "0.35.0/bio/samtools/index"

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
        "-m 4G"
    threads:
        8
    wrapper:
        "0.35.0/bio/samtools/sort"

rule index_final:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        protected("output/mapped/Sample_{sample}.final.bai")
    params:
        "" # optional params string
    wrapper:
        "0.35.0/bio/samtools/index"

rule samtools_flagstat:
    input:
        "output/mapped/Sample_{sample}.final.bam"
    output:
        "output/mapped/Sample_{sample}.final.flagstat"
    wrapper:
        "0.35.0/bio/samtools/flagstat"

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
        rate=$(echo "scale=5;10000000/$nreads" | bc)
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
        "0.35.0/bio/samtools/index"

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
        "macs2 callpeak -t {input} -f BAM -g mm --outdir output/peaks -n {params.sample} 2> {log}"

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
        "macs2 callpeak -t {input} -f BAM -g mm --outdir output/peaks -n {params.sample} 2> {log}"

### Complete QC ###

rule multiqc:
    input:
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
        "0.34.0/bio/multiqc"
