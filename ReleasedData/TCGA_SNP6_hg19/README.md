# CNA profiles of TCGA cases (SNP6) obtained with ASCAT v3

`summary.ascatv3.TCGA.SNP6.penalty70.hg19.tsv` compiles both sample (from TCGA) and genomic (from `ascat.metrics`) information for all high-quality and representative cases.

QC was determined as follows:
- **No ASCAT fit**: ASCAT could not find an optimal ploidy and purity value.
- **Likely normal**: CNA profile matches with a normal sample. Filter: LOH<0.1 & GI<0.1 & purity=1 & ((sex='XX' & ploidy>=1.99 & ploidy<=2.01) | (sex='XY' & ploidy>=1.945 & ploidy<=1.965)).
- **Noisy logR**: logR track (either tumour or normal) is noisy. Filter: tumour_mapd>=0.75 | normal_mapd>=0.75.
- **Likely wrong germline**: matched normal is likely to come from another patient. Filter: tumour_mapd>0.4 & normal_mapd>0.4 & frac_homo>0.1.
- **Oversegmented**: difference in segmentation between logR and BAF tracks. Filter: n_segs_logRBAF_diff>250.
- **Wrong fit**: CNA profile is incorrect, with extreme losses. Filter: mode_majA=0.
- **HD size**: CNA profile contains large segments with homozygous deletion. Filter: homdel_largest>=20e6 | homdel_size>=40e6.
- **Contamination or swap** and **Large CNV**: germline data contains copy-number changes. Filter: we ran aspcf on all germline samples and manually reviewed cases where events were spotted. Cases with large stretches of homozygosity were considered as valid and were not discarded.
- **Pass**: sample has correct metrics.

Only samples defined as **Likely normal** and **Pass** have been included in this data release as the other categories should be discarded from further analyses. Profiles being flagged as **Likely normal** should carefully be interpreted as those do not show evidence of copy-number changes and have a high purity. They could therefore either be true CNA profiles from pure tumours (*e.g.* flat profile for AML) or tumour purity is so low (*e.g.* <20%) that ASCAT fitted the CNA profile on the germline.

After QC, we applied the following strategy to pick one representative sample per case:
- **Pass** were prioritised over **Likely normal** samples.
- Primary tumours were prioritised over recurrences/metastases.
- Blood-derived normals were prioritised over adjacent normal tissues.
- If at least 3 samples were available, we computed Manhattan distances between profiles and selected the two having the minimum distance amongst all samples.
- A consensus quality score was derived based on three features: tumour purity (the higher the better), tumour MAPD (the lower the better) and goodness of fit (the higher the better; GoF). Each feature grants one point and sample with highest score was selected.
- In case of a tie between the two samples, we defined an agreement in MAPD when absolute difference between the two samples was below 0.1 and an agreement in GoF when absolute difference between the two samples was below 5%.
  - If there was an agreement in MAPD but not in GoF, sample with highest GoF was selected.
  - If there was an agreement in GoF but not in MAPD, sample with lowest MAPD was selected.
  - Otherwise, a random sample was selected.