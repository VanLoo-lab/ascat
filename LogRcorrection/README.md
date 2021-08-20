# GC correction file creation

The script *GCfileCreation.sh* can help you create the input for ASCAT’s GC correction for any platform (both array and sequencing data). Please note that this can be quite CPU intensive so should be run on a HPC.

`./GCfileCreation.sh SNPpos SIZES CORES REF`

 - `SNPpos` should contain the list of all probes with their genomic location (tab-separated).

| | Chr | Position |
| --- | --- | --- |
| rs3094315 | 1 | 752566 |
| rs2073813 | 1 | 753541 |
| rs2905040 | 1 | 770216 |
| rs12124819 | 1 | 776546 |
| rs2980314 | 1 | 781258 |

 - `SIZES` contains the length of all reference chromosomes. The annotation for hg19 can be found in the ASCAT repository.
 - `CORES` defines how many cores to use for the computation.
 - `REF` should provide the path to your local reference genome.

The created output file should start like this, with one row per provided probe:

| | Chr | Position | 25bp | 50bp | 100bp | 200bp | 500bp | 1000bp | 2000bp | 5000bp | 10000bp | 20000bp | 50000bp | 100000bp | 200000bp | 500000bp | 1M | 2M | 5M | 10M |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rs3094315 | 1 | 752566 | 0.416667 | 0.44 | 0.36 | 0.36 | 0.406 | 0.435 | 0.428 | 0.4312 | 0.4447 | 0.42725 | 0.43452 | 0.44476 | 0.46216 | 0.483124 | 0.465598 | 0.472959 | 0.51041 | 0.490973 |
| rs2073813 | 1 | 753541 | 0.625 | 0.48 | 0.51 | 0.49 | 0.5 | 0.484 | 0.497 | 0.4414 | 0.4601 | 0.4398 | 0.43728 | 0.4449 | 0.463805 | 0.484478 | 0.465725 | 0.472958 | 0.510401 | 0.490979 |
| rs2905040 | 1 | 770216 | 0.208333 | 0.32 | 0.31 | 0.345 | 0.466 | 0.483 | 0.518 | 0.4914 | 0.4664 | 0.46845 | 0.45344 | 0.44179 | 0.47488 | 0.506966 | 0.471232 | 0.472727 | 0.510679 | 0.49093 |
| rs12124819 | 1 | 776546 | 0.375 | 0.4 | 0.37 | 0.4 | 0.432 | 0.453 | 0.4375 | 0.4488 | 0.4642 | 0.4638 | 0.46552 | 0.44523 | 0.48165 | 0.51042 | 0.475428 | 0.472628 | 0.510835 | 0.490956 |
| rs2980314 | 1 | 781258 | 0.625 | 0.58 | 0.54 | 0.525 | 0.434 | 0.432 | 0.435 | 0.4624 | 0.4711 | 0.46095 | 0.4671 | 0.44497 | 0.48657 | 0.510808 | 0.478539 | 0.472577 | 0.510953 | 0.490999 |

# Replication timing correction file creation
The R script *createReplicTimingfile.R* can help you create the optional replication timing input for ASCAT’s wave correction for any platform (both array and sequencing data). Please note that this currently only works on hg19/GRCh37 due to the lack of publicly available replication data on other genomes. However, one may lift-over hg38/GRCh38 loci to hg19/GRCh37 coordinates to extract information and reset coordinates back to hg38/GRCh38. It also requires that you have the *rtracklayer* package installed (R/Bioconductor).

`./createReplicTimingfile.R SNPpos`

 - `SNPpos` is the same file provided to generate the GC correction file above.

The created output file should start like this, with one row per provided probe:

| | Chr | Position | Bg02es | Bj | Gm06990 | Gm12801 | Gm12812 | Gm12813 | Gm12878 | Helas3 | Hepg2 | Huvec | Imr90 | K562 | Mcf7 | Nhek | Sknsh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rs3094315 | 1 | 752566 | 51.0859069824219 | 59.9739265441895 | 53.6799507141113 | 60.2759857177734 | 52.3710060119629 | 42.9433135986328 | 48.5379791259766 | 53.0320281982422 | 59.2352333068848 | 48.7807922363281 | 67.1048889160156 | 68.2764358520508 | 56.9319763183594 | 55.3978309631348 | 61.0729293823242 |
| rs2073813 | 1 | 753541 | 51.0859069824219 | 59.9739265441895 | 53.6799507141113 | 60.2759857177734 | 52.3710060119629 | 42.9433135986328 | 48.5379791259766 | 53.0320281982422 | 59.2352333068848 | 48.7807922363281 | 67.1048889160156 | 68.2764358520508 | 56.9319763183594 | 55.3978309631348 | 61.0729293823242 |
| rs2905040 | 1 | 770216 | 51.0859069824219 | 59.9739265441895 | 53.6799507141113 | 60.2759857177734 | 52.3710060119629 | 42.9433135986328 | 48.5379791259766 | 53.0320281982422 | 59.2352333068848 | 48.7807922363281 | 67.1048889160156 | 68.2764358520508 | 56.9319763183594 | 55.3978309631348 | 61.0729293823242 |
| rs12124819 | 1 | 776546 | 51.0859069824219 | 59.9739265441895 | 53.6799507141113 | 60.2759857177734 | 52.3710060119629 | 42.9433135986328 | 48.5379791259766 | 53.0320281982422 | 59.2352333068848 | 48.7807922363281 | 67.1048889160156 | 68.2764358520508 | 56.9319763183594 | 55.3978309631348 | 61.0729293823242 |
| rs2980314 | 1 | 781258 | 51.165901184082 | 59.8224296569824 | 53.7252349853516 | 60.2214202880859 | 52.4014015197754 | 42.9043922424316 | 48.5634117126465 | 52.972240447998 | 59.2968635559082 | 48.840877532959 | 67.0139083862305 | 68.229606628418 | 56.9081726074219 | 55.3176765441895 | 61.1264152526855 |
