# Processing WGS data with ASCAT
Although we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) on WGS data to get accurate clonal and subclonal allele-specific copy-number alteration calls, ASCAT can still be used to get a fast ploidy/purity estimate. To this end, we pre-generated a set of loci, alonside with GC content and replication timing correction files.

Briefly, such list was derived from 1000 Genomes Project SNPs ([hg19](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and [hg38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)):

- Biallelic SNPs with allele frequency higher than 0.35 and lower than 0.65 in any population were selected using *BCFtools*.
- Duplicated entries were removed using *R*.
- SNPs located in the [ENCODE blacklisted regions](https://github.com/Boyle-Lab/Blacklist/) were discared.
- SNPs with noisy BAF (distant from 0/0.5/1) in normal samples (a.k.a probloci) as part of the [Battenberg](https://github.com/Wedge-lab/battenberg) package were discarded.

Since hg38 data for the non-PAR region of chrX is not available (as of September 2021), hg38 data for the whole chrX comes from a lift-over from hg19.

GC content and replication timing correction files were then generated using scripts provided in the *[LogRcorrection](../../LogRcorrection)* folder.

Data availability:

- Loci files: [hg19](https://www.dropbox.com/s/l3m0yvyca86lpwb/G1000_loci_hg19.zip) & [hg38](https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip) (unzip and set `alleles.prefix="G1000_loci_hg19_chr"` in `ascat.prepareHTS`)
- Allele files: [hg19](https://www.dropbox.com/s/3fzvir3uqe3073d/G1000_alleles_hg19.zip) & [hg38](https://www.dropbox.com/s/uouszfktzgoqfy7/G1000_alleles_hg38.zip) (unzip and set `loci.prefix="G1000_alleles_hg19_chr"` in `ascat.prepareHTS`)
- GC correction file: [hg19](https://www.dropbox.com/s/v0tgr1esyoh1krw/GC_G1000_hg19.zip) & [hg38](https://www.dropbox.com/s/n7g5dh0ld1hcto8/GC_G1000_hg38.zip) (unzip and set `GCcontentfile="GC_G1000_hg19.txt"` in `ascat.correctLogR`)
- Replication timing correction file: [hg19](https://www.dropbox.com/s/50n7xb06x318tgl/RT_G1000_hg19.zip) & [hg38](https://www.dropbox.com/s/xlp99uneqh6nh6p/RT_G1000_hg38.zip) (unzip and set `replictimingfile="RT_G1000_hg19.txt"` in `ascat.correctLogR`)

Please note that loci files provided above are not 'chr'-based (chromosome names are '1', '2', '3', etc. and not 'chr1', 'chr2', 'chr3', etc.). If your BAMs are 'chr'-based, you will need to add 'chr' (Bash: `for i in {1..22} X; do sed -i 's/^/chr/' G1000_loci_hg19_chr${i}.txt; done`). ASCAT will internally remove 'chr' so the other files (allele, GC correction and RT correction) should not be modified and `chrom_names` (`ascat.prepareHTS`) should be `c(1:22,'X')`.