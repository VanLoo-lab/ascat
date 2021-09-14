# Processing WGS data with ASCAT
Although we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) on WGS data to get accurate clonal and subclonal allele-specific copy-number alteration calls, ASCAT can still be used to get a fast ploidy/purity estimate. To this end, we pre-generated a set of loci, alonside with GC content and replication timing correction files.

Briefly, such list was derived from 1000 Genomes Project SNPs ([hg19](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) and [hg38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)):

- Biallelic SNPs with allele frequency higher than 0.35 and lower than 0.65 in any population were selected using *BCFtools*.
- Duplicated entries were removed using *R*.
- SNPs located in the [ENCODE blacklisted regions](https://github.com/Boyle-Lab/Blacklist/) were discared.
- SNPs with noisy BAF (distant from 0/0.5/1) in normal samples (a.k.a probloci) as part of the [Battenberg](https://github.com/Wedge-lab/battenberg) package were discarded.

Since hg38 data for the non-PAR region of chrX is not available (as of September 2021), hg38 data for the whole chrX comes from a lift-over from hg19.

GC content and replication timing correction files were then generated using scripts provided in the *LogRcorrection* folder.

Data availability:

- Loci files: [hg19]() & [hg38]() (*g1000lociprefix* option in *ascat.prepareHTS*)
- Allele files: [hg19]() & [hg38]() (*g1000allelesprefix* option in *ascat.prepareHTS*)
- GC correction file: [hg19]() & [hg38]() (*GCcontentfile* option in *ascat.GCcorrect*)
- Replication timing correction file: [hg19]() & [hg38]() (*replictimingfile* option in *ascat.GCcorrect*)
