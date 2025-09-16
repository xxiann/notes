---
description: Software for motif discovery and next-gen sequencing analysis
---

# HOMER annotation

Documentation: [http://homer.ucsd.edu/homer/introduction/install.html](http://homer.ucsd.edu/homer/introduction/install.html)

## Installation

```bash
## installation on hpc - atlas8 (need internet access)

mkdir homer
## move the configureHomer.pl into the folder ##

## load required modules
module load gcc
module load zlib
module load samtools

perl configureHomer.pl -install

## add the homer/bin directory to your executable path
export PATH=$PATH:/hpctmp/user/homer/bin/

## to download homer packages
perl /hpctmp/user/homer/configureHomer.pl -list

perl /hpctmp/user/homer/configureHomer.pl -install hg38

## to update 
perl /hpctmp/user/homer/configureHomer.pl -update

# source code
cd /hpctmp/user/anno
export PATH=$PATH:/hpctmp/user/homer/bin/

echo "annotating"

for i in 'all'; do 
    for file in ./${i}/*_merged_peaks.bed; do
        name=$(echo $file|sed 's/_merged_peaks.bed//')
        echo $file
        annotatePeaks.pl \
            $file \
            GRCh38.primary_assembly.genome.fa \
            -gid \
            -gtf gencode.v43.primary_assembly.annotation.gtf \
            -cpu 24 \
            > ${name}.annotatePeaks.txt
    done
done

```



