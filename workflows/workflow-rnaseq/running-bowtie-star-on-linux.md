# running bowtie/STAR on linux

## bowtie

reference: [https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html](https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html)

<pre class="language-sh"><code class="lang-sh">## setting up env
sudo apt install bowtie2
threads=14

<strong>## setting up directories
</strong>mkdir results
mkdir results/counts

## building genome index
echo 'building genome index'
bowtie2-build --threads $threads -f annotation/GRCh38.primary_assembly.genome.fa annotation/GRCh38_index


## aligning
for r1 in $(ls data/*_1.trim.fastq.gz); do
    r2=$(echo $r1|sed 's/_1/_2/')
    prefix=$(echo $r1|sed 's/_1.trim.fastq.gz//'|sed 's/data\///')
    echo $r1 $r2

    echo running bowtie on $prefix
    bowtie2 -p $threads -t \
        -x annotation/GRCh38_index \
        -1 $r1 -2 $r2 \
        -S results/${prefix}.sam

    # sam to bam
    echo 'converting sam to bam'
    samtools view -bS results/${prefix}.sam > results/${prefix}.bam

    rm results/${prefix}.sam

    samtools sort results/${prefix}.bam -o results/${prefix}_sorted.bam

    rm results/${prefix}.bam

    samtools index results/${prefix}_sorted.bam

done

## running feature counts
featureCounts -T $threads -p results/*_sorted.bam \
    -a annotation/gencode.v43.primary_assembly.annotation.gtf \
    -o results/counts/pdx_featurecounts.txt --verbose >featurecount.log


## removing first lline and standardising the name of columns 
for file in $(find ./analysis/counts/ -name '*Counts.txt'); do 
    echo $file
    awk '(NR>1)' $file| \
    sed -e '1s/\(\/mnt\/f\/atacseq_raw\/results\|\/wt\/\|\/rq\/\|\/ko\/\|\/mv411\/\|bowtie2\/merged_library\/\|.mLb.clN.sorted.bam\)//g
' > ${file}.tsv
done
</code></pre>



## STAR

```sh
## setting up directories
mkdir index
mkdir results
mkdir results/counts

## building genome index
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir index/GRCh38_index \
--genomeFastaFiles annotation/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile annotation/gencode.v43.primary_assembly.annotation.gtf \
--sjdbOverhang 74


## aligning
for r1 in $(ls data/*_1.trim.fastq.gz); do
    r2=$(echo $r1|sed 's/_1/_2/')  
    prefix=$(echo $r1|sed 's/_1.trim.fastq.gz//'|sed 's/data\///') 

    STAR --runThreadN 8 \
    --genomeDir index/GRCh38_index \
    --readFilesIn $r1 $r2 \
    --outFileNamePrefix results/$prefix \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI NM MD AS

    ## running feature count; s 0 = unstranded
    featureCounts -T 8 -s 0 results/${prefix}.bam \
    -a annotation/gencode.v43.primary_assembly.annotation.gtf \
    -o results/counts/${name}_featurecounts.txt \

    rm results/${prefix}.bam
 done



# for file in $(ls results/*.bam); do
#     name=$(echo $file|sed 's/.bam//'|sed 's/results\///') 

#     featureCounts -T 8 -s 0 $file \
#     -a annotation/gencode.v43.primary_assembly.annotation.gtf \
#     -o results/counts/${name}_featurecounts.txt \

# done
```
