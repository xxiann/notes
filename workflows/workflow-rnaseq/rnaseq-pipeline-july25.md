---
description: To run on NUS HPC, please change the filename suffix accordingly before use
icon: ballot-check
---

# RNAseq pipeline July'25

You can download the precompiled binary version of subread from [here](https://sourceforge.net/projects/subread/files/subread-2.1.1/) (for 64-bit Linux systems) and upload onto the HPC, unzip it before use `tar zxvf subread-x.x.x.tar.gz`&#x20;

### RNA-seq Pipeline Folder Structure

```markdown
📁 /hpctmp/user/rnaseq
├── 📄 GRCh38_primary_assembly_genome.fa
├── 📄 gencode.v44.primary_assembly.annotation.gtf
├── 📁 rawdata_0/
│ ├── 📁 fastqc_raw/ # FastQC reports for raw reads
│ ├── 📁 trimmed_report/ # Trim Galore reports + temporary trimmed files
│ ├── 📁 trimmed_fq/ # Final trimmed FASTQ files
│ ├── 📁 fastqc_trimmed/ # FastQC reports for trimmed reads
│ ├── 📁 star_alignment/ # STAR intermediate outputs
│ ├── 📁 star_genelevel/ # STAR gene count outputs
│ ├── 📁 bam/ # Final BAM files after STAR alignment
│ ├── 📁 sorted_bam/ # (Unused here – reserved for sorted BAMs if needed)
│ ├── 📄 all_featurecounts.txt # Output from featureCounts
│ ├── 📁 multiqc_report/ # MultiQC summary report
│ └── *.fq.gz # Input FASTQ files (e.g., *_R1.fq.gz, *_R2.fq.gz)

📌 Note:
- The `sorted_bam/` folder is created but not used because STAR already outputs sorted BAMs.
- `star_alignment/` stores temporary STAR outputs before selective files are moved.
- `multiqc_report/` contains the final MultiQC output summarizing the whole run.
```



Done using atlas9

```bash
#!/bin/bash
#PBS -P rnaseq_pipeline
#PBS -q parallel
#PBS -l select=1:ncpus=12:mpiprocs=12:mem=80GB
#PBS -j oe
#PBS -N name

# === SETUP ===
cd "${PBS_O_WORKDIR}"
source /etc/profile.d/rec_modules.sh
source /app1/ebenv

module load STAR/2.7.5b
module load Trim_Galore/0.6.10-GCCcore-11.3.0
module load samtools/1.9

# === VARIABLES ===
RAW_DIR="/hpctmp/user/rnaseq/rawdata_0"
GENOME_DIR="/hpctmp/user/rnaseq"
GENOME_FASTA="${GENOME_DIR}/GRCh38_primary_assembly_genome.fa"
GTF_FILE="${GENOME_DIR}/gencode.v44.primary_assembly.annotation.gtf"

THREADS=12
cd "$RAW_DIR"

# === MAKE OUTPUT DIRS ===
mkdir -p fastqc_raw trimmed_report trimmed_fq fastqc_trimmed star_alignment star_genelevel bam sorted_bam

# === BUILD STAR INDEX (once only) ===
if [ ! -f "${GENOME_DIR}/SAindex" ]; then
  echo "Building STAR genome index..."
  STAR --runThreadN $THREADS \
       --runMode genomeGenerate \
       --genomeDir $GENOME_DIR \
       --genomeFastaFiles $GENOME_FASTA \
       --sjdbGTFfile $GTF_FILE \
       --sjdbOverhang 149
fi

# === MAIN PIPELINE LOOP: replace the ending of the file ===
for R1 in *_R1.fq.gz; do                     #replace the ending of the file
  BASENAME=$(basename "$R1" _R1.fq.gz)       #replace the ending of the file
  R2="${BASENAME}_R2.fq.gz"                  #replace the ending of the file
  # SAMPLE_ID=$(echo "$BASENAME" | cut -d'-' -f2)
  SAMPLE_ID="$BASENAME"

  echo "=== Processing $SAMPLE_ID ==="

  # 1. FastQC on raw reads
  fastqc -o fastqc_raw -t 8 "$R1" "$R2"

  # 2. Trim reads
  trim_galore --paired --fastqc -j 8 --output_dir trimmed_report "$R1" "$R2"

  # 3. Move trimmed FASTQs to their folder
  TRIM_R1="${BASENAME}_R1_val_1.fq.gz"       #replace the ending of the file
  TRIM_R2="${BASENAME}_R2_val_2.fq.gz"       #replace the ending of the file
  mv trimmed_report/"$TRIM_R1" trimmed_fq/
  mv trimmed_report/"$TRIM_R2" trimmed_fq/
  # 4. FastQC on trimmed
  mv trimmed_report/"${TRIM_R1%.fq.gz}_fastqc"* fastqc_trimmed/
  mv trimmed_report/"${TRIM_R2%.fq.gz}_fastqc"* fastqc_trimmed/

  # 5. STAR alignment
  STAR --runThreadN $THREADS \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn trimmed_fq/"$TRIM_R1" trimmed_fq/"$TRIM_R2" \
       --readFilesCommand zcat \
       --outFileNamePrefix star_alignment/"${SAMPLE_ID}_" \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within \
       --quantMode GeneCounts \
       --outSAMattributes NH HI NM MD AS

  # 6. Move output files
  mv star_alignment/"${SAMPLE_ID}_Aligned.out.bam" bam/"${SAMPLE_ID}.bam"
  mv star_alignment/"${SAMPLE_ID}_ReadsPerGene.out.tab" star_genelevel/
  mv star_alignment/"${SAMPLE_ID}_SJ.out.tab" star_genelevel/

  # 7. Sort BAM - already sorted
  # samtools sort -@ 8 -o sorted_bam/"${SAMPLE_ID}.sorted.bam" bam/"${SAMPLE_ID}.bam"
done

# === FEATURECOUNTS ===
export PATH=$PATH:/hpctmp/user/subread-2.1.1-Linux-x86_64/bin/
cd bam
featureCounts -p -T $THREADS -s 0 --countReadPairs \
  -a "$GTF_FILE" \
  -o ../all_featurecounts.txt \
  *.bam

# === MULTIQC ===
module purge
source /app1/ebenv
module load MultiQC/1.14-foss-2022b

cd "$RAW_DIR"
multiqc . -o multiqc_report --interactive -f

echo "=== PIPELINE COMPLETE ==="
```
