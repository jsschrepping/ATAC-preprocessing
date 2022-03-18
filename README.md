Snakemake workflow for pre-processing and peak calling of ATAC-seq data.

Use snakemake/snakemake:v7.2.1 to run the script.

The raw fastq files should be present as one merged fastq.gz file per sample in a fastq/merged_fastq/ directory that is mounted in /data/ in the container.
The bowtie2 index files should be mounted in /data/bowtie2-indexfiles/.