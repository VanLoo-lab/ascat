# GC correction file creation

The script *createGCcontentFile.R* can help you create the input for ASCAT’s GC correction for any platform (both array and sequencing data). Please note that this can be quite CPU/memory intensive so should be run on a HPC. R dependencies: *doParallel*, *foreach* and *Biostrings*.

`Rscript createGCcontentFile.R SNPpos CORES REF`

 - `SNPpos` should contain the list of all probes with their genomic location (tab-separated).

| | Chr | Position |
| --- | --- | --- |
| rs3094315 | 1 | 752566 |
| rs2073813 | 1 | 753541 |
| rs2905040 | 1 | 770216 |
| rs12124819 | 1 | 776546 |
| rs2980314 | 1 | 781258 |

 - `CORES` defines how many cores to use for the computation.
 - `REF` should be the path to your local reference genome. Make sure that chromosome names in `SNPpos` fit with `REF` (e.g. *chr1* is different from *1*). Please note that sequence must be A/C/G/T/N-based (some references contain lowercase letters for specific regions such as repetitive sequences).

The created output file (*GCcontent_SNPloci.txt*) should start like this, with one row per provided probe:

| | Chr | Position | 25bp | 50bp | 100bp | 200bp | 500bp | 1kb | 2kb | 5kb | 10kb | 20kb | 50kb | 100kb | 200kb | 500kb | 1Mb |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rs3094315 | 1 | 752566 | 0.4 | 0.431373 | 0.366337 | 0.363184 | 0.40519 | 0.435564 | 0.428286 | 0.431314 | 0.444656 | 0.427279 | 0.434511 | 0.444756 | 0.462163 | 0.502001 | 0.517332 |
| rs2073813 | 1 | 753541 | 0.6 | 0.490196 | 0.514851 | 0.487562 | 0.499002 | 0.484515 | 0.496752 | 0.441512 | 0.460054 | 0.439778 | 0.437291 | 0.444906 | 0.463808 | 0.50239 | 0.517472 |
| rs2905040 | 1 | 770216 | 0.24 | 0.313725 | 0.306931 | 0.348259 | 0.46507 | 0.483516 | 0.517741 | 0.491502 | 0.466453 | 0.468427 | 0.453431 | 0.441786 | 0.474883 | 0.508137 | 0.522142 |
| rs12124819 | 1 | 776546 | 0.36 | 0.392157 | 0.376238 | 0.39801 | 0.431138 | 0.452547 | 0.437781 | 0.44891 | 0.464254 | 0.463777 | 0.465511 | 0.445236 | 0.481648 | 0.510419 | 0.523123 |
| rs2980314 | 1 | 781258 | 0.6 | 0.568627 | 0.544554 | 0.527363 | 0.433134 | 0.431568 | 0.435282 | 0.462507 | 0.471153 | 0.460927 | 0.467091 | 0.444976 | 0.486573 | 0.510809 | 0.52383 |

# Replication timing correction file creation
The R script *createReplicTimingfile.R* can help you create the optional replication timing input for ASCAT’s wave correction for any platform (both array and sequencing data). R dependency: *rtracklayer*.

`Rscript createReplicTimingFile.R SNPpos`

 - `SNPpos` is the same file provided to generate the GC correction file above. Please note that chromosomes names for replication timing data are *1*, *2*, etc. (not *chr1*, *chr2*, etc.).

The created output file (*ReplicationTiming_SNPloci.txt*) should start like this, with one row per provided probe:

| | Chr | Position | Bg02es | Bj | Gm06990 | Gm12801 | Gm12812 | Gm12813 | Gm12878 | Helas3 | Hepg2 | Huvec | Imr90 | K562 | Mcf7 | Nhek | Sknsh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rs3094315 | 1 | 752566 | 54.905624 | 60.098045 | 58.936962 | 60.200775 | 59.797001 | 57.530655 | 61.163361 | 63.059586 | 68.468796 | 58.805431 | 67.670982 | 67.934761 | 64.411995 | 57.198978 | 58.010719 |
| rs2073813 | 1 | 753541 | 54.966686 | 60.175228 | 58.960648 | 60.191437 | 59.853878 | 57.569115 | 61.162411 | 63.085983 | 68.470589 | 58.821083 | 67.701759 | 67.995796 | 64.448158 | 57.264729 | 57.990158 |
| rs2905040 | 1 | 770216 | 56.000576 | 61.533062 | 59.517155 | 60.315243 | 60.879486 | 58.433681 | 61.421215 | 63.689762 | 68.565933 | 59.307125 | 68.349823 | 69.068848 | 65.135139 | 58.516354 | 57.936165 |
| rs12124819 | 1 | 776546 | 56.495754 | 62.174568 | 59.858677 | 60.523235 | 61.383484 | 58.944229 | 61.69376 | 64.059311 | 68.656227 | 59.664303 | 68.723213 | 69.572342 | 65.509315 | 59.173409 | 58.08889 |
| rs2980314 | 1 | 781258 | 56.79295 | 62.546955 | 60.079067 | 60.68153 | 61.683388 | 59.267509 | 61.892082 | 64.298447 | 68.724228 | 59.908783 | 68.960464 | 69.863716 | 65.745972 | 59.575733 | 58.227413 |

Please note that replication timing information comes from [ENCODE](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/) and is based on hg19. Therefore, script will only works on hg19/GRCh37 SNPs due to the lack of publicly available replication data on other genomes like hg38. However, one may lift-over hg38/GRCh38 loci to hg19/GRCh37 coordinates, extract replication timing information and reset coordinates back to hg38/GRCh38.

# Additional information
Please note that row names in the `SNPpos` file need to be unique but don't need to be SNP IDs (*e.g.* from dbSNP). One can define row names as `${chromosome}_${position}` as long as they are unique and match with row names in logR/BAF files. See our GC correction file for WGS (under [ReferenceFiles/WGS](../ReferenceFiles/WGS)).