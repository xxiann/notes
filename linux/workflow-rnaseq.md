# Workflow for RNA-seq Alignment

<details>
<summary>Structure of folder</summary>
<pre>
|__rnaseq
    |__annotation
        *.fa
        *.gtf
        |__GRCh38_index
    |__data
        *_1.fastq.gz
        *_2.fastq.gz
        |__trimmed
            *_1.trim.fastq.gz
            *_2.trim.fastq.gz
        |__untrimmed
        |__log
    |__fastqc
    |__results
        *.sorted.bam
        *.sorted.bam.bai
        |__counts
            *.txt
            *.logs
            *.tsv
</pre>
</details>


## Environment (Linux - Ubuntu)
```
sudo apt-get update
sudo apt install fastqc
sudo apt install samtools
sudo apt install bowtie2
sudo apt install subread
``` 
## 1. Quality control for fastq
### 1a. Assessing the quality of reads

Running fastqc to generate QC report that con

- format: fastqc -o [output dir] -t [threads] seqfile 1.. seqfileN

```
# start off with your project folder
cd rnaseq
mkdir fastqc

# running fastqc
fastqc -o fastqc -t 8 data/*.fq.gz
```

### 1b. Trimming the reads

There are programs like [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [TrimGalore](https://github.com/FelixKrueger/TrimGalore). The main purpose is to trim low quality reads (usually Q<20) and remove adaptors.

```
# Trimmomatic

TrimmomaticPE -threads 8 data/sample_1.fq.gz data/sample_2.fq.gz \
    data/trimmed/sample_1.trim.fq.gz data/untrimmed/sample_1.untrim.fq.gz \
    data/trimmed/sample_2.trim.fq.gz data/untrimmed/sample_2.untrim.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True SLIDINGWINDOW:4:20 MINLEN:75
```

## 2. Obtaining genome annotations
You can download your preferred genome annotations from the respective genome database. Preferably GRCH38/hg38 version.
- ENSEMBL
- GENCODE
- UCSC
- NCBI 
'explain what .fasta and .gtf are'

```
# start off with your project folder
cd rnaseq
mkdir annotation # create sub foler for annotaion
```
Downloading using url
```
cd annotation 

# Download the fasta file
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz # unmasked file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download the GTF file
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz 
gunzip Homo_sapiens.GRCh38.110.gtf.gz
```

## 3. Gene level quantification

### 3a. Aligning using Bowtie2
```
# start off with your project folder
cd rnaseq
mkdir results 
mkdir results/counts

# set the number of threads used
threads=14
```

Building genome index

- format: bowtie2-build --threads [n] -f [path to annotation file] [prefix to index path]

```
bowtie2-build --threads $threads -f annotation/GRCh38.primary_assembly.genome.fa annotation/GRCh38_index
```

Alignment

```
r1=data/trimmed/sample_1.trim.fastq.gz
r2=data/trimmed/sample_2.trim.fastq.gz

# running bowtie2
bowtie2 -p -t $threads \
    -x annotation/GRCh38_index \
    -1 $r1 -2 $r2 \
    2> data/log/sample_bowtie2.log \
    | samtools view -t $threads -o results/sample.bam

# sorting bam file
samtools sort -@ 6 \
    -o results/sample_sorted.bam
    -T results/sample_sorted
    results/sample.bam

##rm results/sample.bam

# creating bam index file
samtools index results/sample_sorted.bam
```

### 3b. Counting using FeatureCount
```
## running feature count; s 0 = unstranded

featureCounts -T $threads -s 0 results/sample.bam \
    -a annotation/gencode.v43.primary_assembly.annotation.gtf \
    -o results/counts/sample_featurecounts.txt 

#rm results/sample_sorted.bam
```

## 4. Downstream Analysis

SNP? Splice junctions? - requires STAR alignment

Differential Gene Expression Analysis
