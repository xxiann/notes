# nf 2.1.2 script

```bash
#! /bin/bash 	

# "NFCORE_ATACSEQ:ATACSEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES"
samtools faidx GRCh38.primary_assembly.genome.fa
    cut -f 1,2 GRCh38.primary_assembly.genome.fa.fai > GRCh38.primary_assembly.genome.fa.sizes

# NFCORE_ATACSEQ:ATACSEQ:INPUT_CHECK:SAMPLESHEET_CHECK	
check_samplesheet.py \
        samplesheet_wt.csv \
        samplesheet.valid.csv

# NFCORE_ATACSEQ:ATACSEQ:PREPARE_GENOME:GTF2BED	
gtf2bed \
        gencode.v43.primary_assembly.annotation.gtf \
        > gencode.v43.primary_assembly.annotation.bed

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE
[ ! -f  molm13_wt_dmso_REP1_T1_1.fastq.gz ] && ln -s molm13_wt_dmso1_R1.fq.gz molm13_wt_dmso_REP1_T1_1.fastq.gz
        [ ! -f  molm13_wt_dmso_REP1_T1_2.fastq.gz ] && ln -s molm13_wt_dmso1_R2.fq.gz molm13_wt_dmso_REP1_T1_2.fastq.gz
        trim_galore \
            --fastqc \
            --cores 8 \
            --paired \
            --gzip \
            molm13_wt_dmso_REP1_T1_1.fastq.gz \
            molm13_wt_dmso_REP1_T1_2.fastq.gz

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC
printf "%s %s\n" molm13_wt_dmso1_R1.fq.gz molm13_wt_dmso_REP1_T1_1.gz molm13_wt_dmso1_R2.fq.gz molm13_wt_dmso_REP1_T1_2.gz | while read old_name new_name; do
        [ -f "${new_name}" ] || ln -s $old_name $new_name
    done

    fastqc \
        --quiet \
        --threads 6 \
        molm13_wt_dmso_REP1_T1_1.gz molm13_wt_dmso_REP1_T1_2.gz

# NFCORE_ATACSEQ:ATACSEQ:PREPARE_GENOME:GENOME_BLACKLIST_REGIONS
sortBed -i hg38-blacklist.v3.bed -g GRCh38.primary_assembly.genome.fa.sizes | complementBed -i stdin -g GRCh38.primary_assembly.genome.fa.sizes | awk '$1 !~ /chrM/ {print $0}' > GRCh38.include_regions.bed

# NFCORE_ATACSEQ:ATACSEQ:PREPARE_GENOME:GET_AUTOSOMES
get_autosomes.py \
        GRCh38.primary_assembly.genome.fa.fai \
        GRCh38.primary_assembly.genome.fa.autosomes.txt

# NFCORE_ATACSEQ:ATACSEQ:PREPARE_GENOME:TSS_EXTRACT
cat gencode.v43.primary_assembly.annotation.bed | awk -v FS='	' -v OFS='	' '{ if($6=="+") $3=$2+1; else $2=$3-1; print $1, $2, $3, $4, $5, $6;}' > gencode.v43.primary_assembly.annotation.tss.bed

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BOWTIE2_ALIGN	
INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\.rev.1.bt2$//"`
    [ -z "$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\.rev.1.bt2l$//"`
    [ -z "$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \
        -x $INDEX \
        -1 molm13_wt_dmso_REP1_T1_1_val_1.fq.gz -2 molm13_wt_dmso_REP1_T1_2_val_2.fq.gz \
        --threads 12 \
         \
        --rg-id molm13_wt_dmso_REP1_T1 --rg SM:molm13_wt_dmso_REP1 --rg PL:ILLUMINA --rg LB:molm13_wt_dmso_REP1_T1 --rg PU:1 \
        2> molm13_wt_dmso_REP1_T1.Lb.bowtie2.log \
        | samtools view  --threads 12 -o molm13_wt_dmso_REP1_T1.Lb.bam -

    if [ -f molm13_wt_dmso_REP1_T1.Lb.unmapped.fastq.1.gz ]; then
        mv molm13_wt_dmso_REP1_T1.Lb.unmapped.fastq.1.gz molm13_wt_dmso_REP1_T1.Lb.unmapped_1.fastq.gz
    fi

    if [ -f molm13_wt_dmso_REP1_T1.Lb.unmapped.fastq.2.gz ]; then
        mv molm13_wt_dmso_REP1_T1.Lb.unmapped.fastq.2.gz molm13_wt_dmso_REP1_T1.Lb.unmapped_2.fastq.gz
    fi

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT
samtools sort \
         \
        -@ 6 \
        -o molm13_wt_dmso_REP1_T1.Lb.sorted.bam \
        -T molm13_wt_dmso_REP1_T1.Lb.sorted \
        molm13_wt_dmso_REP1_T1.Lb.bam


# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX
samtools \
        index \
        -@ 1 \
         \
        molm13_wt_dmso_REP1_T1.Lb.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:PICARD_MERGESAMFILES_LIBRARY
ln -s molm13_wt_dmso_REP1_T1.Lb.sorted.bam molm13_wt_dmso_REP1.mLb.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_MARKDUPLICATES_PICARD:PICARD_MARKDUPLICATES
picard \
        -Xmx29491M \
        MarkDuplicates \
        --ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
        --INPUT molm13_wt_dmso_REP1.mLb.sorted.bam \
        --OUTPUT molm13_wt_dmso_REP1.mLb.mkD.sorted.bam \
        --REFERENCE_SEQUENCE GRCh38.primary_assembly.genome.fa \
        --METRICS_FILE molm13_wt_dmso_REP1.mLb.mkD.sorted.MarkDuplicates.metrics.txt

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS	
samtools \
        idxstats \
        --threads 0 \
        molm13_wt_dmso_REP1_T1.Lb.sorted.bam \
        > molm13_wt_dmso_REP1_T1.Lb.sorted.bam.idxstats

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS
samtools \
        flagstat \
        --threads 1 \
        molm13_wt_dmso_REP1_T1.Lb.sorted.bam \
        > molm13_wt_dmso_REP1_T1.Lb.sorted.bam.flagstat

# NFCORE_ATACSEQ:ATACSEQ:FASTQ_ALIGN_BOWTIE2:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS	
samtools \
        stats \
        --threads 1 \
        --reference GRCh38.primary_assembly.genome.fa \
        molm13_wt_dmso_REP1_T1.Lb.sorted.bam \
        > molm13_wt_dmso_REP1_T1.Lb.sorted.bam.stats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX
samtools \
        index \
        -@ 1 \
         \
        molm13_wt_dmso_REP1.mLb.mkD.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAMTOOLS_FILTER
samtools view \
        -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 \
        -L GRCh38.include_regions.bed \
        -b molm13_wt_dmso_REP1.mLb.mkD.sorted.bam \
        | bamtools filter \
            -out molm13_wt_dmso_REP1.mLb.flT.sorted.bam \
            -script bamtools_filter_pe.json

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS	
samtools \
        stats \
        --threads 1 \
        --reference GRCh38.primary_assembly.genome.fa \
        molm13_wt_dmso_REP1.mLb.mkD.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.mkD.sorted.bam.stats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT
samtools \
        flagstat \
        --threads 1 \
        molm13_wt_dmso_REP1.mLb.mkD.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.mkD.sorted.bam.flagstat


# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS	
samtools \
        idxstats \
        --threads 0 \
        molm13_wt_dmso_REP1.mLb.mkD.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.mkD.sorted.bam.idxstats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:SAMTOOLS_SORT	
samtools sort \
        -n \
        -@ 6 \
        -o molm13_wt_dmso_REP1.mLb.flT.name_sorted.bam \
        -T molm13_wt_dmso_REP1.mLb.flT.name_sorted \
        molm13_wt_dmso_REP1.mLb.flT.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_REMOVE_ORPHANS
bampe_rm_orphan.py \
            molm13_wt_dmso_REP1.mLb.flT.name_sorted.bam \
            molm13_wt_dmso_REP1.mLb.clN.bam \
            --only_fr_pairs

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT
samtools sort \
         \
        -@ 6 \
        -o molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        -T molm13_wt_dmso_REP1.mLb.clN.sorted \
        molm13_wt_dmso_REP1.mLb.clN.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX
samtools \
        index \
        -@ 1 \
         \
        molm13_wt_dmso_REP1.mLb.clN.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS	
picard \
        -Xmx4915M \
        CollectMultipleMetrics \
        --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
        --INPUT molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        --OUTPUT molm13_wt_dmso_REP1.mLb.clN.CollectMultipleMetrics \
        --REFERENCE_SEQUENCE GRCh38.primary_assembly.genome.fa

# NFCORE_ATACSEQ:ATACSEQ:PICARD_MERGESAMFILES_REPLICATE        	
picard \
            -Xmx29491M \
            MergeSamFiles \
            --SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
            --INPUT molm13_wt_dmso_REP1.mLb.clN.sorted.bam --INPUT molm13_wt_dmso_REP2.mLb.clN.sorted.bam \
            --OUTPUT molm13_wt_dmso.mRp.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT
plotFingerprint \
        --skipZeros --numberOfSamples 500000 --labels molm13_wt_dmso_REP1.mLb.clN \
         \
        --bamfiles molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        --plotFile molm13_wt_dmso_REP1.mLb.clN.plotFingerprint.pdf \
        --outRawCounts molm13_wt_dmso_REP1.mLb.clN.plotFingerprint.raw.txt \
        --outQualityMetrics molm13_wt_dmso_REP1.mLb.clN.plotFingerprint.qcmetrics.txt \
        --numberOfProcessors 12  

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:MACS2_CALLPEAK
macs2 \
        callpeak \
        --keep-dup all --nomodel --pvalue 0.01 \
        --gsize 2861847787.0 \
        --format BAMPE \
        --name molm13_wt_dmso_REP1.mLb.clN \
        --treatment molm13_wt_dmso_REP1.mLb.clN.sorted.bam \

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT	
samtools \
        flagstat \
        --threads 1 \
        molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.clN.sorted.bam.flagstat

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS	
samtools \
        stats \
        --threads 1 \
        --reference GRCh38.primary_assembly.genome.fa \
        molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.clN.sorted.bam.stats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_FILTER_BAM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS
samtools \
        idxstats \
        --threads 0 \
        molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        > molm13_wt_dmso_REP1.mLb.clN.sorted.bam.idxstats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BAM_TO_BIGWIG:BEDTOOLS_GENOMECOV
SCALE_FACTOR=$(grep '[0-9] mapped (' molm13_wt_dmso_REP1.mLb.clN.sorted.bam.flagstat | awk '{print 1000000/$1}')
    echo $SCALE_FACTOR > molm13_wt_dmso_REP1.mLb.clN.scale_factor.txt

    bedtools \
        genomecov \
        -ibam molm13_wt_dmso_REP1.mLb.clN.sorted.bam \
        -bg \
        -scale $SCALE_FACTOR \
        -pc \
         \
    > tmp.bg

    bedtools sort -i tmp.bg > molm13_wt_dmso_REP1.mLb.clN.bedGraph

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:HOMER_ANNOTATEPEAKS	
annotatePeaks.pl \
        molm13_wt_dmso_REP1.mLb.clN_peaks.narrowPeak \
        GRCh38.primary_assembly.genome.fa \
        -gid \
        -gtf gencode.v43.primary_assembly.annotation.gtf \
        -cpu 6 \
        > molm13_wt_dmso_REP1.mLb.clN_peaks.annotatePeaks.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_ATAQV_ATAQV
ataqv \
        --ignore-read-groups \
        --mitochondrial-reference-name chrM \
        --peak-file molm13_wt_dmso_REP1.mLb.clN_peaks.narrowPeak \
        --tss-file gencode.v43.primary_assembly.annotation.tss.bed \
         \
        --autosomal-reference-file GRCh38.primary_assembly.genome.fa.autosomes.txt \
        --metrics-file "molm13_wt_dmso_REP1.ataqv.json" \
        --threads 6 \
        --name molm13_wt_dmso_REP1 \
        NA \
        molm13_wt_dmso_REP1.mLb.mkD.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:FRIP_SCORE
READS_IN_PEAKS=$(intersectBed -a molm13_wt_dmso_REP1.mLb.clN.sorted.bam -b molm13_wt_dmso_REP1.mLb.clN_peaks.narrowPeak -bed -c -f 0.20 | awk -F '	' '{sum += $NF} END {print sum}')
    samtools flagstat molm13_wt_dmso_REP1.mLb.clN.sorted.bam > molm13_wt_dmso_REP1.mLb.clN.sorted.bam.flagstat
    grep 'mapped (' molm13_wt_dmso_REP1.mLb.clN.sorted.bam.flagstat | grep -v "primary" | awk -v a="$READS_IN_PEAKS" -v OFS='	' '{print "molm13_wt_dmso_REP1", a/$1}' > molm13_wt_dmso_REP1.FRiP.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:PLOT_MACS2_QC	
plot_macs2_qc.r \
        -i molm13_wt_va_REP2.mLb.clN_peaks.narrowPeak,molm13_wt_va_REP1.mLb.clN_peaks.narrowPeak,molm13_wt_dmso_REP2.mLb.clN_peaks.narrowPeak,molm13_wt_dmso_REP1.mLb.clN_peaks.narrowPeak \
        -s molm13_wt_va_REP2.mLb.clN,molm13_wt_va_REP1.mLb.clN,molm13_wt_dmso_REP2.mLb.clN,molm13_wt_dmso_REP1.mLb.clN \
        -o ./ -p macs2_peak.mLb.clN

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CONSENSUS_PEAKS:MACS2_CONSENSUS	
sort -T '.' -k1,1 -k2,2n molm13_wt_dmso.mRp.clN_peaks.narrowPeak molm13_wt_va.mRp.clN_peaks.narrowPeak \
        | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > consensus_peaks.mRp.clN.txt

    macs2_merged_expand.py \
        consensus_peaks.mRp.clN.txt \
        molm13_wt_dmso.mRp.clN,molm13_wt_va.mRp.clN \
        consensus_peaks.mRp.clN.boolean.txt \
         \
        --is_narrow_peak

    awk -v FS='	' -v OFS='	' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' consensus_peaks.mRp.clN.boolean.txt > consensus_peaks.mRp.clN.bed

    echo -e "GeneID	Chr	Start	End	Strand" > consensus_peaks.mRp.clN.saf
    awk -v FS='	' -v OFS='	' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' consensus_peaks.mRp.clN.boolean.txt >> consensus_peaks.mRp.clN.saf

    plot_peak_intersect.r -i consensus_peaks.mRp.clN.boolean.intersect.txt -o consensus_peaks.mRp.clN.boolean.intersect.plot.pdf

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CONSENSUS_PEAKS:HOMER_ANNOTATEPEAKS	
annotatePeaks.pl \
        consensus_peaks.mRp.clN.bed \
        GRCh38.primary_assembly.genome.fa \
        -gid \
        -gtf gencode.v43.primary_assembly.annotation.gtf \
        -cpu 6 \
        > consensus_peaks.mRp.clN.annotatePeaks.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CONSENSUS_PEAKS:SUBREAD_FEATURECOUNTS
featureCounts \
        -F SAF -O --fracOverlap 0.2 \
        -p \
        -T 6 \
        -a consensus_peaks.mRp.clN.saf \
        -s 0 \
        -o consensus_peaks.mRp.clN.featureCounts.txt \
        molm13_wt_va_REP2.mLb.clN.sorted.bam molm13_wt_va_REP1.mLb.clN.sorted.bam molm13_wt_dmso_REP2.mLb.clN.sorted.bam molm13_wt_dmso_REP1.mLb.clN.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BAM_TO_BIGWIG:UCSC_BEDGRAPHTOBIGWIG	
bedGraphToBigWig \
        molm13_wt_dmso_REP1.mLb.clN.bedGraph \
        GRCh38.primary_assembly.genome.fa.sizes \
        molm13_wt_dmso_REP1.mLb.clN.bigWig

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS:DEEPTOOLS_COMPUTEMATRIX_SCALE_REGIONS
computeMatrix \
        scale-regions --regionBodyLength 1000 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --missingDataAsZero --skipZeros --smartLabels \
        --regionsFileName gencode.v43.primary_assembly.annotation.bed \
        --scoreFileName molm13_wt_dmso_REP1.mLb.clN.bigWig \
        --outFileName molm13_wt_dmso_REP1.mLb.clN.scale_regions.computeMatrix.mat.gz \
        --outFileNameMatrix molm13_wt_dmso_REP1.mLb.clN.scale_regions.computeMatrix.vals.mat.tab \
        --numberOfProcessors 12

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS:DEEPTOOLS_COMPUTEMATRIX_REFERENCE_POINT	
computeMatrix \
        reference-point --missingDataAsZero --skipZeros --smartLabels --upstream 3000 --downstream 3000 \
        --regionsFileName gencode.v43.primary_assembly.annotation.tss.bed \
        --scoreFileName molm13_wt_dmso_REP1.mLb.clN.bigWig \
        --outFileName molm13_wt_dmso_REP1.mLb.clN.reference_point.computeMatrix.mat.gz \
        --outFileNameMatrix molm13_wt_dmso_REP1.mLb.clN.reference_point.computeMatrix.vals.mat.tab \
        --numberOfProcessors 12

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:MULTIQC_CUSTOM_PEAKS
cat molm13_wt_dmso_REP1.mLb.clN_peaks.narrowPeak | wc -l | awk -v OFS='	' '{ print "molm13_wt_dmso_REP1.mLb.clN_peaks", $1 }' | cat merged_library_peak_count_header.txt - > molm13_wt_dmso_REP1.mLb.clN_peaks.count_mqc.tsv
    cat merged_library_frip_score_header.txt molm13_wt_dmso_REP1.FRiP.txt > molm13_wt_dmso_REP1.mLb.clN_peaks.FRiP_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_MARKDUPLICATES_PICARD:PICARD_MARKDUPLICATES
picard \
        -Xmx29491M \
        MarkDuplicates \
        --ASSUME_SORTED true --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp \
        --INPUT molm13_wt_dmso.mRp.sorted.bam \
        --OUTPUT molm13_wt_dmso.mRp.clN.sorted.bam \
        --REFERENCE_SEQUENCE GRCh38.primary_assembly.genome.fa \
        --METRICS_FILE molm13_wt_dmso.mRp.clN.sorted.MarkDuplicates.metrics.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CONSENSUS_PEAKS:DESEQ2_QC	
deseq2_qc.r \
        --count_file consensus_peaks.mLb.clN.featureCounts.txt \
        --outdir ./ \
        --outprefix consensus_peaks.mLb.clN \
        --cores 6 \
        --id_col 1 --sample_suffix '.mLb.clN.sorted.bam' --count_col 7 --vst TRUE

    sed 's/deseq2_pca/deseq2_pca_1/g' tmp.txt
    sed -i -e 's/DESeq2 /consensus_peaks DESeq2 /g' tmp.txt
    cat tmp.txt consensus_peaks.mLb.clN.pca.vals.txt > consensus_peaks.mLb.clN.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_1/g' tmp.txt
    sed -i -e 's/DESeq2 /consensus_peaks DESeq2 /g' tmp.txt
    cat tmp.txt consensus_peaks.mLb.clN.sample.dists.txt > consensus_peaks.mLb.clN.sample.dists_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_ATAQV_MKARV	
mkarv \
         \
        --concurrency 6 \
        --force \
        ./html/ \
        jsons/*

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CALL_ANNOTATE_PEAKS:PLOT_HOMER_ANNOTATEPEAKS
plot_homer_annotatepeaks.r \
        -i molm13_wt_va_REP2.mLb.clN_peaks.annotatePeaks.txt,molm13_wt_va_REP1.mLb.clN_peaks.annotatePeaks.txt,molm13_wt_dmso_REP2.mLb.clN_peaks.annotatePeaks.txt,molm13_wt_dmso_REP1.mLb.clN_peaks.annotatePeaks.txt \
        -s molm13_wt_va_REP2,molm13_wt_va_REP1,molm13_wt_dmso_REP2,molm13_wt_dmso_REP1 \
        -p macs2_annotatePeaks.mLb.clN \
        -o ./

    find ./ -type f -name "*summary.txt" -exec cat {} \; | cat merged_library_peak_annotation_header.txt - > macs2_annotatePeaks.mLb.clN.summary_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX
samtools \
        index \
        -@ 1 \
         \
        molm13_wt_dmso.mRp.clN.sorted.bam

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:MACS2_CALLPEAK	
macs2 \
        callpeak \
        --keep-dup all --nomodel \
        --gsize 2861847787.0 \
        --format BAMPE \
        --name molm13_wt_dmso.mRp.clN \
        --treatment molm13_wt_dmso.mRp.clN.sorted.bam \

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT 		
samtools \
        flagstat \
        --threads 1 \
        molm13_wt_dmso.mRp.clN.sorted.bam \
        > molm13_wt_dmso.mRp.clN.sorted.bam.flagstat

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS
samtools \
        idxstats \
        --threads 0 \
        molm13_wt_dmso.mRp.clN.sorted.bam \
        > molm13_wt_dmso.mRp.clN.sorted.bam.idxstats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS
samtools \
        stats \
        --threads 1 \
        --reference GRCh38.primary_assembly.genome.fa \
        molm13_wt_dmso.mRp.clN.sorted.bam \
        > molm13_wt_dmso.mRp.clN.sorted.bam.stats

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_BAM_TO_BIGWIG:BEDTOOLS_GENOMECOV
SCALE_FACTOR=$(grep '[0-9] mapped (' molm13_wt_dmso.mRp.clN.sorted.bam.flagstat | awk '{print 1000000/$1}')
    echo $SCALE_FACTOR > molm13_wt_dmso.mRp.clN.scale_factor.txt

    bedtools \
        genomecov \
        -ibam molm13_wt_dmso.mRp.clN.sorted.bam \
        -bg \
        -scale $SCALE_FACTOR \
        -pc \
         \
    > tmp.bg

    bedtools sort -i tmp.bg > molm13_wt_dmso.mRp.clN.bedGraph

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:HOMER_ANNOTATEPEAKS
annotatePeaks.pl \
        molm13_wt_dmso.mRp.clN_peaks.narrowPeak \
        GRCh38.primary_assembly.genome.fa \
        -gid \
        -gtf gencode.v43.primary_assembly.annotation.gtf \
        -cpu 6 \
        > molm13_wt_dmso.mRp.clN_peaks.annotatePeaks.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:FRIP_SCORE
READS_IN_PEAKS=$(intersectBed -a molm13_wt_dmso.mRp.clN.sorted.bam -b molm13_wt_dmso.mRp.clN_peaks.narrowPeak -bed -c -f 0.20 | awk -F '	' '{sum += $NF} END {print sum}')
    samtools flagstat molm13_wt_dmso.mRp.clN.sorted.bam > molm13_wt_dmso.mRp.clN.sorted.bam.flagstat
    grep 'mapped (' molm13_wt_dmso.mRp.clN.sorted.bam.flagstat | grep -v "primary" | awk -v a="$READS_IN_PEAKS" -v OFS='	' '{print "molm13_wt_dmso", a/$1}' > molm13_wt_dmso.FRiP.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:PLOT_MACS2_QC	
plot_macs2_qc.r \
        -i molm13_wt_va.mRp.clN_peaks.narrowPeak,molm13_wt_dmso.mRp.clN_peaks.narrowPeak \
        -s molm13_wt_va.mRp.clN,molm13_wt_dmso.mRp.clN \
        -o ./ -p macs2_peak.mRp.clN

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CONSENSUS_PEAKS:MACS2_CONSENSUS
sort -T '.' -k1,1 -k2,2n molm13_wt_dmso.mRp.clN_peaks.narrowPeak molm13_wt_va.mRp.clN_peaks.narrowPeak \
        | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > consensus_peaks.mRp.clN.txt

    macs2_merged_expand.py \
        consensus_peaks.mRp.clN.txt \
        molm13_wt_dmso.mRp.clN,molm13_wt_va.mRp.clN \
        consensus_peaks.mRp.clN.boolean.txt \
         \
        --is_narrow_peak

    awk -v FS='	' -v OFS='	' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' consensus_peaks.mRp.clN.boolean.txt > consensus_peaks.mRp.clN.bed

    echo -e "GeneID	Chr	Start	End	Strand" > consensus_peaks.mRp.clN.saf
    awk -v FS='	' -v OFS='	' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' consensus_peaks.mRp.clN.boolean.txt >> consensus_peaks.mRp.clN.saf

    plot_peak_intersect.r -i consensus_peaks.mRp.clN.boolean.intersect.txt -o consensus_peaks.mRp.clN.boolean.intersect.plot.pdf

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CONSENSUS_PEAKS:HOMER_ANNOTATEPEAKS	
annotatePeaks.pl \
        consensus_peaks.mLb.clN.bed \
        GRCh38.primary_assembly.genome.fa \
        -gid \
        -gtf gencode.v43.primary_assembly.annotation.gtf \
        -cpu 6 \
        > consensus_peaks.mLb.clN.annotatePeaks.txt

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_CONSENSUS_PEAKS:SUBREAD_FEATURECOUNTS
featureCounts \
        -F SAF -O --fracOverlap 0.2 \
        -p \
        -T 6 \
        -a consensus_peaks.mLb.clN.saf \
        -s 0 \
        -o consensus_peaks.mLb.clN.featureCounts.txt \
        molm13_wt_va_REP2.mLb.clN.sorted.bam molm13_wt_va_REP1.mLb.clN.sorted.bam molm13_wt_dmso_REP2.mLb.clN.sorted.bam molm13_wt_dmso_REP1.mLb.clN.sorted.bam


# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS:DEEPTOOLS_PLOTHEATMAP	
plotHeatmap \
         \
        --matrixFile molm13_wt_dmso_REP1.mLb.clN.reference_point.computeMatrix.mat.gz \
        --outFileName molm13_wt_dmso_REP1.mLb.clN.scale_regions.plotHeatmap.pdf \
        --outFileNameMatrix molm13_wt_dmso_REP1.mLb.clN.scale_regions.plotHeatmap.mat.tab

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_BAM_TO_BIGWIG:UCSC_BEDGRAPHTOBIGWIG	
bedGraphToBigWig \
        molm13_wt_dmso.mRp.clN.bedGraph \
        GRCh38.primary_assembly.genome.fa.sizes \
        molm13_wt_dmso.mRp.clN.bigWig

# NFCORE_ATACSEQ:ATACSEQ:IGV	
find * -type l -name "*.bigWig" -exec echo -e ""{}"\t0,0,178" \; | { grep "^bowtie2/merged_library/bigwig" || test $? = 1; } > mLb_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\t0,0,178" \; | { grep "^bowtie2/merged_library/macs2/narrow_peak" || test $? = 1; } > mLb_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\t0,0,0" \; | { grep "^bowtie2/merged_library/macs2/narrow_peak/consensus" || test $? = 1; } > mLb_bed.igv.txt
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\t0,0,178" \; | { grep "^bowtie2/merged_replicate/bigwig" || test $? = 1; } > mRp_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\t0,0,178" \; | { grep "^bowtie2/merged_replicate/macs2/narrow_peak" || test $? = 1; } > mRp_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\t0,0,0" \; | { grep "^bowtie2/merged_replicate/macs2/narrow_peak/consensus" || test $? = 1; } > mRp_bed.igv.txt

    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/GRCh38.primary_assembly.genome.fa --path_prefix '../../'

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS:DEEPTOOLS_PLOTPROFILE	
plotProfile \
         \
        --matrixFile molm13_wt_dmso_REP1.mLb.clN.scale_regions.computeMatrix.mat.gz \
        --outFileName molm13_wt_dmso_REP1.mLb.clN.reference_point.plotProfile.pdf \
        --outFileNameData molm13_wt_dmso_REP1.mLb.clN.reference_point.plotProfile.tab

# NFCORE_ATACSEQ:ATACSEQ:MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS:DEEPTOOLS_PLOTHEATMAP	
plotHeatmap \
         \
        --matrixFile molm13_wt_dmso_REP2.mLb.clN.reference_point.computeMatrix.mat.gz \
        --outFileName molm13_wt_dmso_REP2.mLb.clN.scale_regions.plotHeatmap.pdf \
        --outFileNameMatrix molm13_wt_dmso_REP2.mLb.clN.scale_regions.plotHeatmap.mat.tab

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CONSENSUS_PEAKS:DESEQ2_QC
deseq2_qc.r \
        --count_file consensus_peaks.mRp.clN.featureCounts.txt \
        --outdir ./ \
        --outprefix consensus_peaks.mRp.clN \
        --cores 6 \
        --id_col 1 --sample_suffix '.mLb.clN.sorted.bam' --count_col 7 --vst TRUE

    sed 's/deseq2_pca/deseq2_pca_1/g' tmp.txt
    sed -i -e 's/DESeq2 /consensus_peaks DESeq2 /g' tmp.txt
    cat tmp.txt consensus_peaks.mRp.clN.pca.vals.txt > consensus_peaks.mRp.clN.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_1/g' tmp.txt
    sed -i -e 's/DESeq2 /consensus_peaks DESeq2 /g' tmp.txt
    cat tmp.txt consensus_peaks.mRp.clN.sample.dists.txt > consensus_peaks.mRp.clN.sample.dists_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:MULTIQC_CUSTOM_PEAKS	
cat molm13_wt_dmso.mRp.clN_peaks.narrowPeak | wc -l | awk -v OFS='	' '{ print "molm13_wt_dmso.mRp.clN_peaks", $1 }' | cat merged_replicate_peak_count_header.txt - > molm13_wt_dmso.mRp.clN_peaks.count_mqc.tsv
    cat merged_replicate_frip_score_header.txt molm13_wt_dmso.FRiP.txt > molm13_wt_dmso.mRp.clN_peaks.FRiP_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MERGED_REPLICATE_CALL_ANNOTATE_PEAKS:PLOT_HOMER_ANNOTATEPEAKS
plot_homer_annotatepeaks.r \
        -i molm13_wt_va.mRp.clN_peaks.annotatePeaks.txt,molm13_wt_dmso.mRp.clN_peaks.annotatePeaks.txt \
        -s molm13_wt_va,molm13_wt_dmso \
        -p macs2_annotatePeaks.mRp.clN \
        -o ./

    find ./ -type f -name "*summary.txt" -exec cat {} \; | cat merged_replicate_peak_annotation_header.txt - > macs2_annotatePeaks.mRp.clN.summary_mqc.tsv

# NFCORE_ATACSEQ:ATACSEQ:MULTIQC
multiqc \
        -f \
        --title "atacseq" \
         \
        .

######################################################################################################################
```
