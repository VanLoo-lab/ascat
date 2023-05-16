# Processing WES data with ASCAT

The following files have been derived from the reference files in the [Battenberg package](https://github.com/Wedge-lab/battenberg). No filter related to the genomic location was applied so they contain both exonic, intronic and intergenic SNPs. When using such files, one **must** either 1) provide a BED file which defines regions of interest (`BED_file` in `ascat.prepareHTS`) or 2) downsample the files so they are tailored to your sequencing design (option 2 will speed-up ASCAT but only recommended for advanced users). This is because a WES experiment would only cover a small fraction of the genome so we have to provide an exhaustive list of SNPs to start with, considering that only a subset would be covered. As such, reference files for WGS have a lower resolution (there are not meant to be downsampled) and must not be used for processing WES data. Also, because reference files for WES contain an exhaustive list of SNPs, they must not be used for processing WGS. 

Please note that such files can also be used for processing targeted sequencing data (with an appropriate BED file). Since they contain exonic/intronic/intergenic SNPs, they should be applicable to a broad range of designs.

Data availability:

- Loci files: [hg19](https://drive.google.com/file/d/1R2MON6M77kh3M2v7nSOM1_MOEXXnpIbH/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1BqqYUFdQ3uVllBwSkBYs-jrRJFsmB-No/view?usp=share_link) (unzip and set `alleles.prefix="G1000_loci_hg19_chr"` in `ascat.prepareHTS`)
- Allele files: [hg19](https://drive.google.com/file/d/1NtwlPiArWeFaMO1PjuH7aNhhNZ3bMX5J/view?usp=share_link) & [hg38](https://drive.google.com/file/d/101tHEt4r81b-bZxCE50edRXuOCh4bo-M/view?usp=share_link) (unzip and set `loci.prefix="G1000_alleles_hg19_chr"` in `ascat.prepareHTS`)
- GC correction file: [hg19](https://drive.google.com/file/d/1-WZGH0XTtWFjyJV6ucWQ-yaoYtHcztvq/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1-0wOc5woK27tEpUpSNsz0U3xwI-xvBb3/view?usp=share_link) (unzip and set `GCcontentfile="GC_G1000_hg19.txt"` in `ascat.correctLogR`)
- Replication timing correction file: [hg19](https://drive.google.com/file/d/1unad1TkrXV-EFcRM8KBUCXpLXWCIxWyv/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1Etl5HYdC4TNw0oj4AxvbJR_1a1atYBFt/view?usp=share_link) (unzip and set `replictimingfile="RT_G1000_hg19.txt"` in `ascat.correctLogR`)

### File format

#### Loci file
One file per chromosome, no header, the first column is the chromosome name and the second column is the position.

| | |
| --- | --- |
| 1 | 10642 |
| 1 | 11008 |
| 1 | 11012 |

#### Allele file
One file per chromosome, one header, the first column is the position and the second and third columns are reference and alternate nucleotides with the following conversion: A=1, C=2, G=3 and T=4.

| position | a0 | a1 |
| --- | --- | --- |
| 10642 | 3 | 1 |
| 11008 | 2 | 3 |
| 11012 | 2 | 3 |

In the example above, SNP at position 10642 is G>A and SNPs at 11008 and 11012 are both C>G.

#### GC correction file
One single file, one header, the first column is the SNP ID, the second/third columns are chromosome/position and the other columns are GC% around SNPs with different window sizes.

| | Chr | Position | 25bp | 50bp | 100bp | 200bp | 500bp | 1kb | 2kb | 5kb | 10kb | 20kb | 50kb | 100kb | 200kb | 500kb | 1Mb |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1_10642 | 1 | 10642 | 0.92 | 0.823529 | 0.762376 | 0.761194 | 0.722555 | 0.677323 | 0.625457 | 0.595799 | 0.590039 | 0.5845710.533734 | 0.458927 | 0.421891 | 0.425195 | 0.423964 |
| 1_11008 | 1 | 11008 | 0.72 | 0.745098 | 0.722772 | 0.741294 | 0.730539 | 0.705295 | 0.594703 | 0.593501 | 0.594541 | 0.5832120.534297 | 0.457987 | 0.42164 | 0.425088 | 0.423964 |
| 1_11012 | 1 | 11012 | 0.76 | 0.705882 | 0.742574 | 0.741294 | 0.726547 | 0.706294 | 0.595202 | 0.593964 | 0.594478 | 0.5831820.53433 | 0.457971 | 0.421633 | 0.425084 | 0.423964 |

#### Replication timing file
One single file, one header, the first column is the SNP ID, the second/third columns are chromosome/position and the other columns are replication timing data in different cell lines.

| | Chr | Position | Bg02es | Bj | Gm06990 | Gm12801 | Gm12812 | Gm12813 | Gm12878 | Helas3 | Hepg2 | Huvec | Imr90 | K562 | Mcf7 | Nhek | Sknsh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1_10642 | 1 | 10642 | 49.509453 | 62.858498 | 52.757858 | 61.294971 | 51.757736 | 43.72905 | 48.088467 | 54.11837 | 58.062084 | 47.565636 | 68.790581 | 68.970825 | 57.467934 | 56.897934 | 60.012413 |
| 1_11008 | 1 | 11008 | 49.509453 | 62.858498 | 52.757858 | 61.294971 | 51.757736 | 43.72905 | 48.088467 | 54.11837 | 58.062084 | 47.565636 | 68.790581 | 68.970825 | 57.467934 | 56.897934 | 60.012413 |
| 1_11012 | 1 | 11012 | 49.509453 | 62.858498 | 52.757858 | 61.294971 | 51.757736 | 43.72905 | 48.088467 | 54.11837 | 58.062084 | 47.565636 | 68.790581 | 68.970825 | 57.467934 | 56.897934 | 60.012413 |

### 'chr'-based versus non 'chr'-based reference

Please note that loci files provided above are not 'chr'-based (chromosome names are '1', '2', '3', etc. and not 'chr1', 'chr2', 'chr3', etc.). If your BAMs are 'chr'-based, you will need to add 'chr' (Bash: `for i in {1..22} X; do sed -i 's/^/chr/' G1000_loci_hg19_chr${i}.txt; done`). ASCAT will internally remove 'chr' so the other files (allele, GC correction and RT correction) should not be modified and `chrom_names` (`ascat.prepareHTS`) should be `c(1:22,'X')`.