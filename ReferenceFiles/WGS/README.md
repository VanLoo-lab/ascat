# Processing WGS data with ASCAT
Although we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) on WGS data to get accurate clonal and subclonal allele-specific copy-number alteration calls, ASCAT can still be used to get a fast ploidy/purity estimate. To this end, we pre-generated a set of loci, alongside with GC content and replication timing correction files.

Briefly, such a list was derived from 1000 Genomes Project SNPs ([hg19](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and [hg38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)):

- Biallelic SNPs with allele frequency higher than 0.35 and lower than 0.65 in any population were selected using *BCFtools*.
- Duplicated entries were removed using *R*.
- SNPs located in the [ENCODE blacklisted regions](https://github.com/Boyle-Lab/Blacklist/) were discarded.
- SNPs with noisy BAF (distant from 0/0.5/1) in normal samples (a.k.a probloci) as part of the [Battenberg](https://github.com/Wedge-lab/battenberg) package were discarded (*probloci_270415.txt.gz* for hg19 and *probloci.zip* for hg38).

Since hg38 data for the non-PAR region of chrX is not available (as of September 2021), hg38 data for the whole chrX comes from a lift-over from hg19.

GC content and replication timing correction files were then generated using scripts provided in the *[LogRcorrection](../../LogRcorrection)* folder.

Data availability:

- Loci files: [hg19](https://drive.google.com/file/d/1W0Yxkj9osFs6raEc18TuNhBG9wthnPoB/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1uD2emA-LRJKYabrKSOH-XL40QnyJW1Lh/view?usp=share_link) (unzip and set `alleles.prefix="G1000_loci_hg19_chr"` in `ascat.prepareHTS`)
- Allele files: [hg19](https://drive.google.com/file/d/1ztA_LrBVsyMJJ6niiqcf_5uvTvxJPOze/view?usp=share_link) & [hg38](https://drive.google.com/file/d/14iDvfUegZ4eF5wSuxkeEA7yR-yZYkLjt/view?usp=share_link) (unzip and set `loci.prefix="G1000_alleles_hg19_chr"` in `ascat.prepareHTS`)
- GC correction file: [hg19](https://drive.google.com/file/d/1JB4tBGJmzmjYDtpYY9estRxdTYQSLldD/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1919sBMW5_ul8dXHNgy58_dFy79fpfGY8/view?usp=share_link) (unzip and set `GCcontentfile="GC_G1000_hg19.txt"` in `ascat.correctLogR`)
- Replication timing correction file: [hg19](https://drive.google.com/file/d/1K1qSS8NUzMM8sXdVD8wan9wQ0WRoW593/view?usp=share_link) & [hg38](https://drive.google.com/file/d/1coDFcPp3bCWr-9ZaNoCzQUwewdEfH7fg/view?usp=share_link) (unzip and set `replictimingfile="RT_G1000_hg19.txt"` in `ascat.correctLogR`)

Please note that loci files provided above are not 'chr'-based (chromosome names are '1', '2', '3', etc. and not 'chr1', 'chr2', 'chr3', etc.). If your BAMs are 'chr'-based, you will need to add 'chr' (Bash: `for i in {1..22} X; do sed -i 's/^/chr/' G1000_loci_hg19_chr${i}.txt; done`). ASCAT will internally remove 'chr' so the other files (allele, GC correction and RT correction) should not be modified and `chrom_names` (`ascat.prepareHTS`) should be `c(1:22,'X')`.