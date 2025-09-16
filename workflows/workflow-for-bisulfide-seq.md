# workflow for bisulfide-seq

genome annotation

* [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/)

Bismark

* [https://github.com/FelixKrueger/Bismark/releases/tag/v0.24.2](https://github.com/FelixKrueger/Bismark/releases/tag/v0.24.2)
* [https://felixkrueger.github.io/Bismark/bismark/genome\_preparation/](https://felixkrueger.github.io/Bismark/bismark/genome_preparation/)

## quality control

[https://felixkrueger.github.io/Bismark/bismark/library\_types/](https://felixkrueger.github.io/Bismark/bismark/library_types/)

* needing to trim 5’ and 3’ for illumina aligned + deduplicate after alignment

using trim\_galore [https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim\_Galore\_User\_Guide.md](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

## alignment

*   non\_directional: The sequencing library was constructed in a non strand-specific manner, alignments to all four bisulfite strands will be reported. Default: OFF.

    * (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary to the original strands are merely theoretical and should not exist in reality. Specifying directional alignments (which is the default) will only run 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel and report these alignments. This is the recommended option for sprand-specific libraries).

    ```
    Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
    Setting parallelization to single-threaded (default)

    Summary of all aligner options: -q -N 1 --score-min L,0,-0.2 -p 8 --reorder --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500
    ```

    ```
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /home/svu/user/project/bisulfide/non_dir/
    Setting parallelization to single-threaded (default)

    Summary of all aligner options: -q -N 1 --score-min L,0,-0.2 -p 4 --reorder --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500
    ```
* in this case the alignment should be directional, since none was aligned to the complimentary strand

## deduplication

* [https://felixkrueger.github.io/Bismark/bismark/library\_types/](https://felixkrueger.github.io/Bismark/bismark/library_types/)
* to remove PCR bias caused by overamplification

## methylation extraction

* [https://felixkrueger.github.io/Bismark/bismark/methylation\_extraction/](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/)

## differential methylation analysis ([tutorial](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#id16))

* using [**methylKit**](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#3_Comparative_analysis)
* uses Bismark `.bismark.cov.gz` as input
* can be done on 3 different levels
  * differential methylated bases
  * differential methylated regions (CpG islands, or specified regions)
  * differential methylated binned regions
* differential methylated bases
  1. loading data → set design → filter by coverage → normalise by coverage size → merge data (extract consensus bases) → further filtering based on SD (if more than 2 samples) → (pool to 2 groups: control vs treatment - for Fisher test) → differential analysis → filter and obtain significant one → annotate → export as bed files with %methylation\_difference info for visualisation on UCSC
