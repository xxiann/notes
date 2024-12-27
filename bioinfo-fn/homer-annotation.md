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
export PATH=$PATH:/hpctmp/xiaoxian/homer/bin/

## to download homer packages
perl /hpctmp/xiaoxian/homer/configureHomer.pl -list

perl /hpctmp/xiaoxian/homer/configureHomer.pl -install hg38

## to update 
perl /hpctmp/xiaoxian/homer/configureHomer.pl -update
```



