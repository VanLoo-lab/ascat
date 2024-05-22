#' Function to plot logR and BAF
#' @noRd
plotLogRandBAF=function(Plotdir, SNPpos, LogR, BAF) {
  myLegend=function(ch) {
    abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = ch[[j]]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos, 2, chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5, lty=1, col="lightgrey")
    }
  }
  stopifnot(identical(rownames(LogR), rownames(BAF)))
  stopifnot(identical(rownames(LogR), rownames(SNPpos)))
  stopifnot(identical(colnames(LogR), colnames(BAF)))
  samples=colnames(LogR)
  chrs=as.vector(unique(SNPpos[, 1]))
  last=0
  ch=list()
  SNPorder=vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke=SNPpos[SNPpos[, 1]==chrs[i], ]
    chrpos=chrke[, 2]
    names(chrpos)=rownames(chrke)
    chrpos=sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))
    SNPorder[ch[[i]]]=names(chrpos)
    last=last+length(chrpos)
  }
  SNPpos=SNPpos[SNPorder, ]
  LogR=LogR[SNPorder, , drop=FALSE]
  BAF=BAF[SNPorder, , drop=FALSE]
  for (i in 1:dim(LogR)[2]) {
    png(filename=paste(Plotdir, "/", samples[i], ".png", sep=""), width=2000, height=1000, res=200)
    par(mar=c(0.5, 5, 5, 0.5), mfrow=c(2, 1), cex=0.4, cex.main=3, cex.axis=2, pch=20)
    # plot logR
    plot(c(1, dim(LogR)[1]), c(-2, 2), type="n", xaxt="n", main=paste0(samples[i], ", LogR"), xlab="", ylab="")
    points(LogR[, i], col="red")
    points(LogR[, i], col="#77000011")
    myLegend(ch)
    # plot mirrored BAF
    plot(c(1, dim(BAF)[1]), c(0, 1), type="n", xaxt="n", main=paste0(samples[i], ", BAF"), xlab="", ylab="")
    plot_as_minor_or_major <- runif(length(BAF[, i]))<0.5
    points(ifelse(plot_as_minor_or_major, BAF[, i], 1-BAF[, i]), col="red")
    points(ifelse(plot_as_minor_or_major, BAF[, i], 1-BAF[, i]), col="#77000011")
    myLegend(ch)
    dev.off()
  }
}

#' Function to compute logR based on counts (normals only)
#' @noRd
computeLogR=function(NORMAL_COUNTS_tot) {
  logR=data.frame(apply(NORMAL_COUNTS_tot, 2, function(x) {
    x=log2(x)
    x=x-mean(x)
    return(x)
  }), check.names=FALSE)
  logRref = apply(logR, 1, median)
  logR=data.frame(apply(logR, 2, function(x) {
    x = x-logRref
    x = x-mean(x)
    return(x)
  }), check.names=FALSE)
  return(logR)
}

#' Generate a cleaned list of SNPs from normal samples.
#'
#' @param Worksheet A table with the following columns: Patient_ID, Normal_ID, Normal_file and Gender.
#' @param Workdir The folder where output should go.
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr").
#' @param minCounts Minimum depth, in normal samples, required for a SNP to be considered.
#' @param is_chr_based A boolean indicating whether data is "chr"-based (e.g. 'chr1' instead of '1'). Default=FALSE.
#' @param X_nonPAR Vector containing genomic coordinates (start & stop) of nonPAR region on X. Default=NULL.
#' @param chrom_names A vector containing the names of chromosomes to be considered (optional, default=c(1:22, "X")).
#' @param plotQC A boolean to generate QC reports as PNGs (optional, default=TRUE).
#' @noRd
getLociFromNormals=function(Worksheet, Workdir, alleles.prefix, minCounts, is_chr_based=FALSE, X_nonPAR=NULL, chrom_names=c(1:22, "X"), plotQC=TRUE) {
  stopifnot((is.null(X_nonPAR)) || (length(X_nonPAR)==2 && all(is.numeric(X_nonPAR))))
  min_samples_nonPAR=10 # This defines how many females we require to filter the nonPAR region
  # Read all alleleCount files (counts>=0)
  print("      Getting allele counts...")
  NORMAL_COUNTS=foreach(INDEX=1:nrow(Worksheet)) %dopar% {
    INDEX=get("INDEX")
    return(readAlleleCountFiles(paste0(Workdir, "/alleleCounts/", Worksheet$Patient_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "_unfiltered_chr"), ".txt", chrom_names, 0, keep_chr_string=is_chr_based))
  }
  names(NORMAL_COUNTS)=Worksheet$Normal_ID
  stopifnot(all(sapply(2:length(NORMAL_COUNTS), function(x) identical(rownames(NORMAL_COUNTS[[1]]), rownames(NORMAL_COUNTS[[x]])))))
  print("      Getting allelic information...")
  allele_data=readAllelesFiles(alleles.prefix, ".txt", chrom_names, add_chr_string=is_chr_based)
  stopifnot(identical(rownames(NORMAL_COUNTS[[1]]), rownames(allele_data)))
  # Get ref/alt/tot information for all cases based on alleleCounts and allele data
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) {
    x=x[, 3:6]
    x$ref=x[cbind(1:nrow(allele_data), allele_data$a0)]
    x$alt=x[cbind(1:nrow(allele_data), allele_data$a1)]
    x$tot=x$ref+x$alt
    x$vaf=x$alt/x$tot
    x=x[, c("ref", "alt", "tot", "vaf")]
    return(x)
  })
  # Flag samples having too few SNPs (counts>minCounts)
  COVERED=sapply(NORMAL_COUNTS, function(x) length(which(x$tot>=minCounts)))
  if (plotQC) {
    png(paste0(Workdir, "/plotQC/Low_number_of_SNPs.png"), height=15, width=15, units="cm", res=300, pointsize=6)
    par(mar=c(8.4, 4.2, 1.05, 1.05))
    barplot(sort(COVERED), names.arg=names(sort(COVERED)), space=0, las=3, cex.names=0.5, ylab=paste0("Number of SNPs (counts>=", minCounts, ")"), cex.lab=1.25)
    abline(h=median(COVERED)/2, col="red")
    dev.off()
  }
  TO_REMOVE=names(which(COVERED<median(COVERED)/2))
  if (length(TO_REMOVE)>0) {
    print(paste0("   Remove samples with very low number of SNPs above threshold: ", paste(TO_REMOVE, collapse=", ")))
    NORMAL_COUNTS=NORMAL_COUNTS[setdiff(names(NORMAL_COUNTS), TO_REMOVE)]
    Worksheet=Worksheet[-which(Worksheet$Normal_ID %in% TO_REMOVE), ]
  }
  rm(TO_REMOVE, COVERED)
  # Flag samples having low coverages
  COVERED=sapply(NORMAL_COUNTS, function(x) sum(x$tot[x$tot>=minCounts]))
  if (plotQC) {
    png(paste0(Workdir, "/plotQC/Low_coverage.png"), height=15, width=15, units="cm", res=300, pointsize=6)
    par(mar=c(8.4, 4.2, 1.05, 1.05))
    barplot(sort(COVERED), names.arg=names(sort(COVERED)), space=0, las=3, cex.names=0.5, ylab=paste0("Total coverage (counts>=", minCounts, ")"), cex.lab=1.25)
    abline(h=median(COVERED)/2, col="red")
    dev.off()
  }
  TO_REMOVE=names(which(COVERED<median(COVERED)/2))
  if (length(TO_REMOVE)>0) {
    print(paste0("   Remove samples with very low coverage: ", paste(TO_REMOVE, collapse=", ")))
    NORMAL_COUNTS=NORMAL_COUNTS[setdiff(names(NORMAL_COUNTS), TO_REMOVE)]
    Worksheet=Worksheet[-which(Worksheet$Normal_ID %in% TO_REMOVE), ]
  }
  rm(TO_REMOVE, COVERED)
  # Get SNPs covered in all samples
  IDs=rownames(allele_data)[apply(do.call(cbind, lapply(NORMAL_COUNTS, function(x) x$tot)), 1, min)>0]
  print(paste0("   Keep SNPs covered in all samples: ", length(IDs)))
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[IDs, ])
  allele_data=allele_data[IDs, ]
  rm(IDs)
  #######################################################################################
  # FILTER: keep SNPs with enough coverage (>=minCounts in at least 90% of cases) #
  #######################################################################################
  TO_KEEP=rowSums(do.call(cbind, lapply(NORMAL_COUNTS, function(x) x$tot))>=minCounts)>=length(NORMAL_COUNTS)*0.9
  print(paste0("   Keep SNPs with enough coverage in most cases: ", length(which(TO_KEEP)), " (-", round((nrow(allele_data)-length(which(TO_KEEP)))/nrow(allele_data)*100, 2), "%)"))
  allele_data=allele_data[TO_KEEP, ]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[TO_KEEP, ])
  rm(TO_KEEP)
  ############################################
  # FILTER: keep SNPs with finite BAF values #
  ############################################
  TO_KEEP=apply(do.call(cbind, lapply(NORMAL_COUNTS, function(x) x[, "vaf", drop=FALSE])), 1, function(x) all(is.finite(x)))
  print(paste0("   Keep SNPs with finite BAF values: ", length(which(TO_KEEP)), " (-", round((nrow(allele_data)-length(which(TO_KEEP)))/nrow(allele_data)*100, 2), "%)"))
  allele_data=allele_data[TO_KEEP, ]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[TO_KEEP, ])
  rm(TO_KEEP)
  #################################################
  # FILTER: remove homozygous SNPs in all samples #
  #################################################
  # The probabilistic method (based on binom.test) can be very slow because it needs to be computed for nSNPs x nSamples.
  # We can speed-up the process by:
  # 1) generating the unique list of alt/tot counts (so if several SNPs have 10/50, then it's computed once)
  # 2) processing the 0-0.5 BAF space (two SNPs with 10/50 and 40/50 will have similar confidence intervals as conf and 1-rev(conf))
  all_alt_tot=unique(do.call(rbind, lapply(NORMAL_COUNTS, function(x) x[, c("alt", "tot")])))
  all_alt_tot_inf=all_alt_tot[which(all_alt_tot$alt<=all_alt_tot$tot/2), ]
  all_alt_tot_sup=all_alt_tot[which(all_alt_tot$alt>all_alt_tot$tot/2), ]
  all_alt_tot_sup$alt=all_alt_tot_sup$tot-all_alt_tot_sup$alt
  all_alt_tot=rbind(all_alt_tot_inf, all_alt_tot_sup)
  rm(all_alt_tot_inf, all_alt_tot_sup)
  all_alt_tot=unique(all_alt_tot)
  all_alt_tot=all_alt_tot[order(all_alt_tot$tot, all_alt_tot$alt), ]
  rownames(all_alt_tot)=paste0(all_alt_tot$alt, "_", all_alt_tot$tot)
  # Get confidence intervals
  conf=foreach(MAX=unique(all_alt_tot$tot), .combine=c) %dopar% {
    MAX=get("MAX")
    lapply(all_alt_tot$alt[all_alt_tot$tot==MAX], function(x) {
      binom.test(x=x, n=MAX, alternative = "two.sided", p=0.01, conf.level = 0.98)$conf.int
    })
  }
  # Get the interpretation
  all_alt_tot$interpretation=sapply(conf, function(x) {
    if (x[1]<=0.5 && x[2]>=0.5 && (x[1]<=0.01 || x[2]>=0.99)) return("noisy") # This will occur for low-coverage SNPs with mixed signal
    if (x[1]<=0.5 && x[2]>=0.5) return("het")
    if (x[1]<=0.01 || x[2]>=0.99) return("hom")
    return("noisy")
  })
  # Now, assign genotypes to SNPs in all samples
  GENOTYPES=foreach(INDEX=1:length(NORMAL_COUNTS), .combine=cbind) %dopar% {
    INDEX=get("INDEX")
    all_alt_tot[paste0(ifelse(NORMAL_COUNTS[[INDEX]]$alt>NORMAL_COUNTS[[INDEX]]$tot/2, NORMAL_COUNTS[[INDEX]]$tot-NORMAL_COUNTS[[INDEX]]$alt, NORMAL_COUNTS[[INDEX]]$alt), "_", NORMAL_COUNTS[[INDEX]]$tot), "interpretation"]
  }
  rm(conf)
  TO_REMOVE=which(apply(GENOTYPES, 1, function(x) all(x=="hom")))
  if (length(which(Worksheet$Gender=="XX"))<min_samples_nonPAR) {
    # Here, we don't have enough females for the nonPAR region. We rescue SNPs on the nonPAR region if they are not noisy in males
    index_nonPAR=which(allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), "X") & allele_data$position>=X_nonPAR[1] & allele_data$position<=X_nonPAR[2])
    index_nonPAR=index_nonPAR[which(apply(GENOTYPES[index_nonPAR, which(Worksheet$Gender=="XY")], 1, function(x) all(x=="hom")))]
    TO_REMOVE=setdiff(TO_REMOVE, index_nonPAR)
    rm(index_nonPAR)
  }
  if (length(TO_REMOVE)>0) {
    print(paste0("   Remove homozygous SNPs in all samples: ", nrow(allele_data)-length(TO_REMOVE), " (-", round(length(TO_REMOVE)/nrow(allele_data)*100, 2), "%)"))
    allele_data=allele_data[-TO_REMOVE, ]
    NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[-TO_REMOVE, ])
    GENOTYPES=GENOTYPES[-TO_REMOVE, ]
  }
  rm(all_alt_tot, TO_REMOVE)
  #################################################################
  # FILTER: remove SNPs with noisy BAF values in multiple samples #
  #################################################################
  NORMAL_COUNTS_tot=do.call(cbind, lapply(NORMAL_COUNTS, function(x) x$tot))
  INDEX_autosomes=which(allele_data$chromosome %in% paste0(ifelse(is_chr_based, "chr", ""), 1:22) | (allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), "X") & (allele_data$position<X_nonPAR[1] | allele_data$position>X_nonPAR[2])))
  TO_KEEP_autosomes=which(rowSums(NORMAL_COUNTS_tot[INDEX_autosomes, ]>=minCounts & (GENOTYPES[INDEX_autosomes, ]=="noisy"))<2.5*0.053876*rowSums(NORMAL_COUNTS_tot[INDEX_autosomes, ]>=minCounts&GENOTYPES[INDEX_autosomes, ] %in% c("noisy", "het")))
  TO_KEEP_autosomes=INDEX_autosomes[TO_KEEP_autosomes]
  INDEX_XX=which(Worksheet$Gender=="XX")
  INDEX_XY=which(Worksheet$Gender=="XY")
  INDEX_nonPAR=which(allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), "X") & allele_data$position>=X_nonPAR[1] & allele_data$position<=X_nonPAR[2])
  if (length(INDEX_nonPAR)>0) {
    if (length(INDEX_XX)>=min_samples_nonPAR) {
      # Here, we do have enough female samples to proceed so we want SNPs to be clean (low noisy/(noisy+het) ratio)
      TO_KEEP_nonPAR=which(rowSums(NORMAL_COUNTS_tot[INDEX_nonPAR, INDEX_XX]>=minCounts & (GENOTYPES[INDEX_nonPAR, INDEX_XX]=="noisy"))<2.5*0.053876*rowSums(NORMAL_COUNTS_tot[INDEX_nonPAR, INDEX_XX]>=minCounts&GENOTYPES[INDEX_nonPAR, INDEX_XX] %in% c("noisy", "het")))
    } else {
      # Here, there are too few females so consider all SNPs as possible candidates
      TO_KEEP_nonPAR=1:length(INDEX_nonPAR)
    }
    # SNPs must be clean in males (all hom)
    TO_KEEP_nonPAR=TO_KEEP_nonPAR[apply(GENOTYPES[INDEX_nonPAR[TO_KEEP_nonPAR], INDEX_XY], 1, function(x) all(x=="hom"))]
    TO_KEEP_nonPAR=INDEX_nonPAR[TO_KEEP_nonPAR]
  } else {
    TO_KEEP_nonPAR=NA
  }
  TO_KEEP=sort(unique(c(TO_KEEP_autosomes, TO_KEEP_nonPAR)))
  print(paste0("   Remove SNPs with noisy BAF: ", length(TO_KEEP), " (-", round((nrow(allele_data)-length(TO_KEEP))/nrow(allele_data)*100, 2), "%)"))
  allele_data=allele_data[TO_KEEP, ]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[TO_KEEP, ])
  GENOTYPES=GENOTYPES[TO_KEEP, ]
  rm(TO_KEEP, NORMAL_COUNTS_tot, INDEX_autosomes, INDEX_nonPAR, INDEX_XX, INDEX_XY, TO_KEEP_autosomes, TO_KEEP_nonPAR)
  ##################################################################
  # FILTER: remove SNPs close to each other with similar genotypes #
  ##################################################################
  GENOTYPES[GENOTYPES=="noisy"]="het"
  SNP_dist=diff(allele_data[, 2])
  TO_REMOVE = sapply(which(SNP_dist<=75), function(poske) {
    hetp1 = ifelse(GENOTYPES[poske, ]=="het", 1, 0)
    hetp2 = ifelse(GENOTYPES[poske+1, ]=="het", 1, 0)
    # if less than 10% of the calls are different
    if (sum(abs(hetp1-hetp2))<0.05 * (sum(hetp1)+sum(hetp2))) {
      if (SNP_dist[max(poske-1, 1)]<SNP_dist[min(poske+1, length(SNP_dist))]) {
        return(poske)
      } else {
        return(poske+1)
      }
    } else {
      return(NA)
    }
  })
  TO_REMOVE = unique(na.omit(TO_REMOVE))
  if (length(TO_REMOVE)>0) {
    print(paste0("   Keep SNPs separated (distance and genotype): ", nrow(allele_data)-length(TO_REMOVE), " (-", round(length(TO_REMOVE)/nrow(allele_data)*100, 2), "%)"))
    allele_data=allele_data[-TO_REMOVE, ]
    NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[-TO_REMOVE, ])
  }
  rm(TO_REMOVE, SNP_dist, GENOTYPES)
  #######################################
  # FILTER: remove SNPs with noisy logR #
  #######################################
  # Compute logR
  TOT=do.call(cbind, lapply(NORMAL_COUNTS, function(x) x$tot))
  INDEX_XY=which(Worksheet$Gender=="XY")
  INDEX_nonPAR=which(allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), "X") & allele_data$position>=X_nonPAR[1] & allele_data$position<=X_nonPAR[2])
  if (length(INDEX_XY)>0 && length(INDEX_nonPAR)>0) TOT[INDEX_nonPAR, INDEX_XY]=TOT[INDEX_nonPAR, INDEX_XY]*2
  logR=computeLogR(TOT)
  rm(TOT, INDEX_nonPAR, INDEX_XY)
  # Get standard deviation
  sdkes = apply(logR, 1, sd)
  TO_KEEP=sdkes<=2*median(sdkes)
  if (plotQC) {
    myRLE=rle(allele_data$chromosome)
    png(paste0(Workdir, "/plotQC/LogR_noise.png"), height=7.5, width=15, units="cm", res=300, pointsize=6)
    par(mar=c(2.1, 4.2, 1.05, 1.05), lwd=0.5)
    plot(NULL, xlim=c(1, length(sdkes)), ylim=c(0, max(c(max(sdkes), 2*median(sdkes)))), ylab="LogR SD", xlab="", xaxt="n", xaxs="i")
    abline(h=2*median(sdkes), col="red", lwd=0.5)
    abline(v=cumsum(myRLE$lengths)[-length(myRLE$lengths)], lwd=0.5)
    points(sdkes, pch=16, cex=0.125, col=ifelse(TO_KEEP, "black", "red"))
    axis(1, at=cumsum(c(0, myRLE$lengths[-length(myRLE$lengths)]))+myRLE$lengths/2, labels=myRLE$values, las=2, cex.axis=0.5)
    dev.off()
    rm(myRLE)
  }
  print(paste0("   Remove SNPs with noisy logR: ", length(which(TO_KEEP)), " (-", round((nrow(allele_data)-length(which(TO_KEEP)))/nrow(allele_data)*100, 2), "%)"))
  allele_data=allele_data[TO_KEEP, ]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS, function(x) x[TO_KEEP, ])
  rm(TO_KEEP, sdkes, logR)
  # Re-compute cleaned logR
  TOT=do.call(cbind, lapply(NORMAL_COUNTS, function(x) x[, "tot", drop=FALSE]))
  INDEX_XY=which(Worksheet$Gender=="XY")
  INDEX_nonPAR=which(allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), "X") & allele_data$position>=X_nonPAR[1] & allele_data$position<=X_nonPAR[2])
  if (length(INDEX_XY)>0 && length(INDEX_nonPAR)>0) TOT[INDEX_nonPAR, INDEX_XY]=TOT[INDEX_nonPAR, INDEX_XY]*2
  logR=computeLogR(TOT)
  rm(TOT, INDEX_nonPAR, INDEX_XY)
  colnames(logR)=names(NORMAL_COUNTS)
  # Flag noisy sample(s)
  sdsam = apply(logR, 2, sd)
  if (plotQC) {
    png(paste0(Workdir, "/plotQC/Noisy_logR.png"), height=15, width=15, units="cm", res=300, pointsize=6)
    par(mar=c(8.4, 4.2, 1.05, 1.05))
    barplot(sort(sdsam), names.arg=names(sort(sdsam)), space=0, las=3, cex.names=0.5, ylab="Noise in logR", cex.lab=1.25, ylim=c(0, max(1.5*median(sdsam), max(sdsam))))
    abline(h=1.5*median(sdsam), col="red")
    dev.off()
  }
  IDs=names(which(sdsam>1.5*median(sdsam)))
  if (length(IDs)>0) {
    warning(paste0("Noisy sample", if (length(IDs)>1) "s", ": ", paste0(IDs, collapse=", ")))
  }
  rm(IDs, sdsam)
  # Write & plot logR+BAF values
  if (!dir.exists(paste0(Workdir, "/Normal_data"))) dir.create(paste0(Workdir, "/Normal_data"))
  write.table(cbind(allele_data[, 1:2], do.call(cbind, lapply(NORMAL_COUNTS, function(x) x$vaf))), paste0(Workdir, "/Normal_data/Normal_BAF.txt"), sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
  write.table(cbind(allele_data[, 1:2], logR), paste0(Workdir, "/Normal_data/Normal_LogR.txt"), sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
  BAF=do.call(cbind, lapply(NORMAL_COUNTS, function(x) x[, "vaf", drop=FALSE]))
  colnames(BAF)=names(NORMAL_COUNTS)
  plotLogRandBAF(paste0(Workdir, "/Normal_data"), allele_data[, 1:2], logR, BAF)
  # Write output
  if (!dir.exists(paste0(Workdir, "/alleleData/Cleaned"))) dir.create(paste0(Workdir, "/alleleData/Cleaned"))
  for (CHR in chrom_names) {
    tmp=allele_data[which(allele_data[, 1]==paste0(ifelse(is_chr_based, "chr", ""), CHR)), ]
    write.table(tmp[, 1:2], paste0(Workdir, "/alleleData/Cleaned/loci_chr", CHR, ".txt"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    write.table(tmp[, 2:4], paste0(Workdir, "/alleleData/Cleaned/alleleData_chr", CHR, ".txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  }; rm(CHR, tmp)
}

#' Method to extract a curated list of SNPs covered by a targeted sequencing experiment.
#'
#' From a complete set of loci (alleles.prefix), this method will keep SNPs falling into the targeted design (based on BED_file) and check allele counts in normal samples (listed in Worksheet). The cleaned list of loci/allele files will be located under Workdir/alleleData/Cleaned/.
#'
#' @param Worksheet A tab-separated file with the following columns: Patient_ID, Normal_ID, Normal_file and Gender (additional columns can be provided but will not be used). Must contain one single normal per patient. Normal_file can either be paths to BAMs/CRAMs or paths to pre-computed (zipped) alleleCounts (e.g. "sample_alleleCounts_chr"). Gender must either be XX (females) or XY (males).
#' @param Workdir The folder where output should go (will be created if it doesn't exist).
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr").
#' @param BED_file A BED file for only looking at SNPs within specific intervals. Must fit with the design used for targeted sequencing.
#' @param allelecounter_exe Path to the allele counter executable.
#' @param genomeVersion Genome version, either 'hg19' or 'hg38'.
#' @param nthreads The number of parallel processes to speed up the process (optional, default=1).
#' @param minCounts Minimum depth required in the normal for a SNP to be considered (optional, default=10).
#' @param is_chr_based A boolean indicating whether data is "chr"-based (e.g. 'chr1' instead of '1'; optional, default=FALSE).
#' @param chrom_names A vector containing the names of chromosomes to be considered (optional, default=c(1:22, "X")). Do not set it to paste0("chr", c(1:22, "X")) if data is "chr"-based.
#' @param min_base_qual Minimum base quality required for a read to be counted (optional, default=20).
#' @param min_map_qual Minimum mapping quality required for a read to be counted (optional, default=35).
#' @param ref.fasta FASTA file used for generating CRAMs (optional, default=NA).
#' @param plotQC A boolean to generate QC reports as PNGs (optional, default=TRUE).
#' @export
ascat.prepareTargetedSeq=function(Worksheet, Workdir, alleles.prefix, BED_file, allelecounter_exe, genomeVersion, nthreads=1,
                                  minCounts=10, is_chr_based=FALSE, chrom_names=c(1:22, "X"), min_base_qual=20, min_map_qual=35, ref.fasta=NA, plotQC=TRUE) {
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  requireNamespace("foreach")
  requireNamespace("doParallel")

  ########
  # Init #
  ########
  stopifnot(genomeVersion %in% c("hg19", "hg38"))
  if (genomeVersion=="hg19") {
    X_nonPAR=c(2699521, 154931043)
  } else if (genomeVersion=="hg38") {
    X_nonPAR=c(2781480, 155701382)
  }
  Worksheet=read.table(Worksheet, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  stopifnot(nrow(Worksheet)>1)
  stopifnot(all(c("Patient_ID", "Normal_ID", "Normal_file", "Gender") %in% colnames(Worksheet)))
  stopifnot(length(which(duplicated(Worksheet$Patient_ID)))==0)
  stopifnot(length(which(duplicated(Worksheet$Sample_ID)))==0)
  stopifnot(all(Worksheet$Gender %in% c("XX", "XY")))
  registerDoParallel(cores=nthreads)
  if (!dir.exists(Workdir)) dir.create(Workdir)
  if (!dir.exists(paste0(Workdir, "/alleleData/Raw"))) dir.create(paste0(Workdir, "/alleleData/Raw"), recursive=TRUE)
  if (!dir.exists(paste0(Workdir, "/alleleCounts"))) dir.create(paste0(Workdir, "/alleleCounts"))
  if (!dir.exists(paste0(Workdir, "/plotQC")) && plotQC) dir.create(paste0(Workdir, "/plotQC"))
  Process_HTS_file=FALSE
  Process_AC_file=FALSE

  # Two scenarios here: either Normal_file are BAMs/CRAMs or paths to pre-computed (zipped) alleleCounts
  if (all(sapply(strsplit(Worksheet$Normal_file, "\\."), function(x) x[length(x)]) %in% c("bam", "BAM", "cram", "CRAM"))) {
    stopifnot(all(file.exists(Worksheet$Normal_file)))
    Process_HTS_file=TRUE
  } else if (all(sapply(Worksheet$Normal_file, function(x) all(file.exists(paste0(x, chrom_names, ".txt")))))) {
    Process_AC_file=TRUE
    SUFFIX=".txt"
  } else if (all(sapply(Worksheet$Normal_file, function(x) all(file.exists(paste0(x, chrom_names, ".txt.gz")))))) {
    Process_AC_file=TRUE
    SUFFIX=".txt.gz"
  } else {
    stop("Worksheet$Normal_file seems incorrect; Please check that input files exist and have the right extension (bam/BAM/cram/CRAM from HTS data or txt/txt.gz for pre-computed allele counts.)")
  }

  # subset SNPs based on BED
  print("Subsetting SNPs based on BED")
  allele_data=readAllelesFiles(alleles.prefix, ".txt", chrom_names, add_chr_string=is_chr_based)
  BED=read.table(BED_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1:3]
  colnames(BED)=c("chr", "start", "end")
  BED$start=BED$start+1 # Start is 0-based in BED files
  BED=BED[BED$chr %in% paste0(ifelse(is_chr_based, "chr", ""), chrom_names), ]
  if (nrow(BED)==0) stop("Major issue with BED file, please double-check its content")
  overlaps=findOverlaps(GRanges(seqnames=BED$chr, ranges=IRanges(start=BED$start, end=BED$end)),
                        GRanges(seqnames=allele_data$chromosome, ranges=IRanges(start=allele_data$position, end=allele_data$position)))
  allele_data=allele_data[unique(overlaps@to), ]
  for (chr in chrom_names) {
    tmp=allele_data[allele_data$chromosome==paste0(ifelse(is_chr_based, "chr", ""), chr), ]
    tmp=tmp[order(tmp$position), ]
    write.table(tmp[, -1], file=paste0(Workdir, "/alleleData/Raw/alleleData_chr", chr, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(tmp[, 1:2], file=paste0(Workdir, "/alleleData/Raw/loci_chr", chr, ".txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    rm(tmp)
  }; rm(chr, allele_data, BED, overlaps)

  if (Process_AC_file) {
    allele_data=readAllelesFiles(paste0(Workdir, "/alleleData/Raw/alleleData_chr"), ".txt", chrom_names, add_chr_string=is_chr_based)
    allele_data=lapply(paste0(ifelse(is_chr_based, "chr", ""), chrom_names), function(x) rownames(allele_data[allele_data$chromosome==x, ]))
    names(allele_data)=paste0(ifelse(is_chr_based, "chr", ""), chrom_names)
  }

  for (INDEX in 1:nrow(Worksheet)) {
    print(paste0("   Processing normal sample: ", Worksheet$Normal_ID[INDEX], " (", Worksheet$Patient_ID[INDEX], "; ", INDEX, "/", nrow(Worksheet), ")"))
    if (!dir.exists(paste0(Workdir, "/alleleCounts/", Worksheet$Patient_ID[INDEX], "/", Worksheet$Normal_ID[INDEX]))) dir.create(paste0(Workdir, "/alleleCounts/", Worksheet$Patient_ID[INDEX], "/", Worksheet$Normal_ID[INDEX]), recursive=TRUE)
    if (Process_HTS_file) {
      # Run alleleCounter on all samples (normals only)
      foreach(CHR=chrom_names) %dopar% {
        CHR=get("CHR")
        ascat.getAlleleCounts(seq.file=Worksheet$Normal_file[INDEX],
                              output.file=paste0(Workdir, "/alleleCounts/", Worksheet$Patient_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "_unfiltered_chr", CHR, ".txt"),
                              loci.file=paste0(Workdir, "/alleleData/Raw/loci_chr", CHR, ".txt"),
                              min.base.qual=min_base_qual,
                              min.map.qual=min_map_qual,
                              allelecounter.exe=allelecounter_exe,
                              ref.fasta=ref.fasta)
      }
    } else {
      foreach(CHR=chrom_names) %dopar% {
        CHR=get("CHR")
        LOCI=fread(paste0(Worksheet$Normal_file[INDEX], CHR, SUFFIX), sep="\t", data.table=FALSE, showProgress=FALSE, quote="")
        rownames(LOCI)=paste0(LOCI[, 1], "_", LOCI[, 2])
        stopifnot(all(allele_data[[paste0(ifelse(is_chr_based, "chr", ""), CHR)]] %in% rownames(LOCI)))
        LOCI=LOCI[allele_data[[paste0(ifelse(is_chr_based, "chr", ""), CHR)]], ]
        write.table(LOCI, file=paste0(Workdir, "/alleleCounts/", Worksheet$Patient_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "/", Worksheet$Normal_ID[INDEX], "_unfiltered_chr", CHR, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      }
    }
  }; rm(INDEX)

  ###########################
  # 3) Clean normal samples #
  ###########################
  print("Clean normal samples")
  getLociFromNormals(Worksheet=Worksheet,
                     Workdir=Workdir,
                     alleles.prefix=paste0(Workdir, "/alleleData/Raw/alleleData_chr"),
                     chrom_names=chrom_names,
                     minCounts=minCounts,
                     is_chr_based=is_chr_based,
                     X_nonPAR=X_nonPAR,
                     plotQC=plotQC)
}
