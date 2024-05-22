#' Obtain allele counts for a given set of loci through external program alleleCounter.
#'
#' @param seq.file A BAM/CRAM alignment file on which the counter should be run.
#' @param output.file The file where output should go.
#' @param loci.file A file with SNP loci.
#' @param min.base.qual The minimum base quality required for it to be counted (optional, default=20).
#' @param min.map.qual The minimum mapping quality required for it to be counted (optional, default=35).
#' @param allelecounter.exe A pointer to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @param ref.fasta A FASTA file for CRAM processing (optional).
#' @author sd11, tl
#' @export
ascat.getAlleleCounts = function(seq.file, output.file, loci.file, min.base.qual=20, min.map.qual=35, allelecounter.exe="alleleCounter", ref.fasta=NA) {
  if (!file.exists(seq.file) || file.info(seq.file)$size==0) {warning("seq.file does not seem to exist or is empty"); return()}
  if (!file.exists(loci.file) || file.info(loci.file)$size==0) {warning("loci.file does not seem to exist or is empty"); return()}
  cmd = paste(allelecounter.exe,
              "-b", seq.file,
              "-l", loci.file,
              "-o", output.file,
              "-m", min.base.qual,
              "-q", min.map.qual)
  # alleleCounter >= v4.0.0 is sped up considerably when run in dense-snp mode
  counter_version = system(paste(allelecounter.exe, "--version"), intern = TRUE)
  if (as.integer(substr(x = counter_version, start = 1, stop = 1)) >= 4) {
    cmd = paste(cmd, "--dense-snps")
  }
  # if ref.fasta is provided, add it to the command line
  if (!is.na(ref.fasta)) {
    cmd = paste(cmd, "-r", ref.fasta)
  }
  EXIT_CODE=system(cmd, wait=TRUE)
  stopifnot(EXIT_CODE==0)
}

#' Obtain BAF and LogR from the allele counts.
#'
#' @param samplename String, name of the sample.
#' @param tumourAlleleCountsFile.prefix Prefix of the allele counts files for the tumour (e.g. "Tumour_alleleFrequencies_chr").
#' @param normalAlleleCountsFile.prefix Prefix of the allele counts files for the normal (e.g. "Normal_alleleFrequencies_chr").
#' @param tumourLogR_file File where LogR from the tumour will be written.
#' @param tumourBAF_file File where BAF from the tumour will be written.
#' @param normalLogR_file File where LogR from the normal will be written.
#' @param normalBAF_file File where BAF from the normal will be written.
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr")
#' @param gender Gender information, either 'XX' (=female) or 'XY' (=male).
#' @param genomeVersion Genome version, either 'hg19' or 'hg38'.
#' @param chrom_names A vector with allowed chromosome names (optional, default=c(1:22, 'X')). Do not set it to paste0('chr', c(1:22, 'X')) if data is 'chr'-based.
#' @param minCounts Minimum depth, in normal samples, required for a SNP to be considered (optional, default=20).
#' @param BED_file A BED file for only looking at SNPs within specific intervals (optional, default=NA).
#' @param probloci_file A file (chromosome <tab> position; no header) containing specific loci to ignore (optional, default=NA).
#' @param seed A seed to be set when randomising the alleles (optional, default=as.integer(Sys.time())).
#' @author dw9, sd11, tl
#' @export
ascat.getBAFsAndLogRs = function(samplename, tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file, alleles.prefix, gender, genomeVersion, chrom_names=c(1:22, "X"), minCounts=20, BED_file=NA, probloci_file=NA, seed=as.integer(Sys.time())) {
  set.seed(seed)
  stopifnot(gender %in% c("XX", "XY"))
  stopifnot(genomeVersion %in% c("hg19", "hg38"))
  # Load data, only keep SNPs with enough coverage
  tumour_input_data = readAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chrom_names, 1)
  normal_input_data = readAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chrom_names, minCounts)
  allele_data = readAllelesFiles(alleles.prefix, ".txt", chrom_names)
  # Synchronise DFs
  matched_data = Reduce(intersect, list(rownames(tumour_input_data), rownames(normal_input_data), rownames(allele_data)))
  tumour_input_data = tumour_input_data[rownames(tumour_input_data) %in% matched_data, ]
  normal_input_data = normal_input_data[rownames(normal_input_data) %in% matched_data, ]
  allele_data = allele_data[rownames(allele_data) %in% matched_data, ]
  rm(matched_data)
  # If a probloci file is provided, remove those
  if (!is.na(probloci_file)) {
    stopifnot(file.exists(probloci_file) && file.info(probloci_file)$size>0)
    probloci=data.frame(fread(probloci_file, sep="\t", showProgress=FALSE, header=TRUE), stringsAsFactors=FALSE)
    probloci=paste0(gsub("^chr", "", probloci[, 1]), "_", probloci[, 2])
    probloci=which(rownames(tumour_input_data) %in% probloci)
    if (length(probloci)>0) {
      tumour_input_data = tumour_input_data[-probloci, ]
      normal_input_data = normal_input_data[-probloci, ]
      allele_data = allele_data[-probloci, ]
    } else {
      warning("The probloci did not remove any SNPs, it might be worth checking the data.")
    }
    rm(probloci)
  }
  stopifnot(isTRUE(all.equal(allele_data[, 1], tumour_input_data[, 1])) && isTRUE(all.equal(allele_data[, 1], normal_input_data[, 1])))
  stopifnot(isTRUE(all.equal(allele_data[, 2], tumour_input_data[, 2])) && isTRUE(all.equal(allele_data[, 2], normal_input_data[, 2])))
  tumour_input_data = tumour_input_data[, 3:6]
  normal_input_data = normal_input_data[, 3:6]
  # If a BED is provided, only look at SNPs within those intervals
  if (!is.na(BED_file)) {
    stopifnot(file.exists(BED_file) && file.info(BED_file)$size>0)
    BED=read.table(BED_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1:3]
    colnames(BED)=c("chr", "start", "end")
    BED$chr=gsub("^chr", "", BED$chr)
    BED$start=BED$start+1 # Start is 0-based in BED files
    BED=BED[BED$chr %in% chrom_names, ]
    if (nrow(BED)==0) stop("Major issue with BED file, please double-check its content")
    requireNamespace("GenomicRanges")
    requireNamespace("IRanges")
    overlaps=findOverlaps(GRanges(seqnames=BED$chr, ranges=IRanges(start=BED$start, end=BED$end)),
                          GRanges(seqnames=allele_data$chromosome, ranges=IRanges(start=allele_data$position, end=allele_data$position)))
    if (length(overlaps)>0) {
      tumour_input_data=tumour_input_data[unique(overlaps@to), ]
      normal_input_data=normal_input_data[unique(overlaps@to), ]
      allele_data=allele_data[unique(overlaps@to), ]
    } else {
      print(head(allele_data))
      print(head(BED))
      stop("The overlap between the BED file and loci is empty. Data must be checked!")
    }
    rm(BED, overlaps)
  }
  # Obtain depth for both alleles for tumour and normal
  len = nrow(allele_data)
  tumour_input_data$REF=tumour_input_data[cbind(1:len, allele_data[, 3])]
  tumour_input_data$ALT=tumour_input_data[cbind(1:len, allele_data[, 4])]
  normal_input_data$REF=normal_input_data[cbind(1:len, allele_data[, 3])]
  normal_input_data$ALT=normal_input_data[cbind(1:len, allele_data[, 4])]
  # Make sure that ALT+REF fit with minimal counts
  TO_KEEP=which(tumour_input_data$REF+tumour_input_data$ALT>=1 & normal_input_data$REF+normal_input_data$ALT>=minCounts)
  stopifnot(length(TO_KEEP)>0)
  allele_data=allele_data[TO_KEEP, ]
  tumour_input_data=tumour_input_data[TO_KEEP, ]
  normal_input_data=normal_input_data[TO_KEEP, ]
  rm(TO_KEEP)
  # Prepare allele counts to derive BAF and logR
  len = nrow(allele_data)
  mutCount1 = tumour_input_data$REF
  mutCount2 = tumour_input_data$ALT
  totalTumour = mutCount1 + mutCount2
  normCount1 = normal_input_data$REF
  normCount2 = normal_input_data$ALT
  totalNormal = normCount1 + normCount2
  rm(tumour_input_data, normal_input_data)
  normalBAF = vector(length=len, mode="numeric")
  tumourBAF = vector(length=len, mode="numeric")
  normalLogR = vector(length=len, mode="numeric")
  tumourLogR = vector(length=len, mode="numeric")
  # Output raw (=unmirrored) BAF from some downstream analyses (e.g. refphase)
  normalBAF_unmirrored=normCount2/totalNormal
  tumourBAF_unmirrored=mutCount2/totalTumour
  germline.BAF_unmirrored = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=normalBAF_unmirrored, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  tumor.BAF_unmirrored = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=tumourBAF_unmirrored, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  colnames(tumor.BAF_unmirrored)[3]=samplename
  colnames(germline.BAF_unmirrored)[3]=samplename
  write.table(tumor.BAF_unmirrored, file=gsub("\\.txt$", "_rawBAF.txt", tumourBAF_file), row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
  write.table(germline.BAF_unmirrored, file=gsub("\\.txt$", "_rawBAF.txt", normalBAF_file), row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
  rm(normalBAF_unmirrored, tumourBAF_unmirrored, germline.BAF_unmirrored, tumor.BAF_unmirrored)
  # Randomise A and B alleles
  selector = round(runif(len))
  normalBAF[which(selector==0)] = normCount1[which(selector==0)] / totalNormal[which(selector==0)]
  normalBAF[which(selector==1)] = normCount2[which(selector==1)] / totalNormal[which(selector==1)]
  tumourBAF[which(selector==0)] = mutCount1[which(selector==0)] / totalTumour[which(selector==0)]
  tumourBAF[which(selector==1)] = mutCount2[which(selector==1)] / totalTumour[which(selector==1)]
  # Normalise tumourLogR to normalLogR
  tumourLogR = totalTumour/totalNormal
  tumourLogR = log2(tumourLogR/mean(tumourLogR, na.rm=TRUE))
  rm(selector)
  # For males, chrX needs to be adjusted as logR baseline will be 0 because of T/N ratio
  if (gender=="XY") {
    # PAR1 and PAR2 information should be a mix of chrX and chrY so we should expect 1+1 (1 from X and 1 from Y).
    # nonPAR should be X-specific and baseline is 1+0 so logR needs to be decreased according to gamma parameter (ascat.runAscat)
    if (genomeVersion=="hg19") {
      nonPAR=c(2699521, 154931043)
    } else if (genomeVersion=="hg38") {
      nonPAR=c(2781480, 155701382)
    }
    nonPAR=which(allele_data$chromosome %in% c("X", "chrX") & allele_data$position>=nonPAR[1] & allele_data$position<=nonPAR[2])
    tumourLogR[nonPAR]=tumourLogR[nonPAR]-1
  }
  # Create the output data.frames
  tumor.LogR = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, logr=tumourLogR, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  tumor.BAF = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=tumourBAF, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  germline.LogR = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, logr=normalLogR, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  germline.BAF = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=normalBAF, ID=rownames(allele_data), row.names=4, stringsAsFactors=FALSE)
  colnames(tumor.LogR)[3]=samplename
  colnames(tumor.BAF)[3]=samplename
  colnames(germline.LogR)[3]=samplename
  colnames(germline.BAF)[3]=samplename
  # Save data.frames to disk
  write.table(tumor.LogR, file=tumourLogR_file, row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
  write.table(tumor.BAF, file=tumourBAF_file, row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
  write.table(germline.LogR, file=normalLogR_file, row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
  write.table(germline.BAF, file=normalBAF_file, row.names=TRUE, quote=FALSE, sep="\t", col.names=NA)
}

#' Synchronise SNPs across files
#'
#' @param samplename String, name of the sample.
#' @param tumourLogR_file File where LogR from the tumour will be read and overwritten.
#' @param tumourBAF_file File where BAF from the tumour will be read and overwritten.
#' @param normalLogR_file File where LogR from the normal will be read and overwritten.
#' @param normalBAF_file File where BAF from the normal will be read and overwritten.
#' @author tl
#' @export
ascat.synchroniseFiles=function(samplename, tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file) {
  # read all files
  FILES=lapply(c(tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file), function(x) {
    tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=TRUE, na.strings=c("-Inf", "Inf", "NA", "NaN", "", "-")), row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    colnames(tmp)=c("Chromosome", "Position", samplename)
    tmp=tmp[!is.na(tmp[, 3]), ]
    return(tmp)
  })
  names(FILES)=c(tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file)
  # get IDs shared between DFs
  IDs=Reduce(intersect, lapply(FILES, rownames))
  FILES=lapply(FILES, function(x) x[rownames(x) %in% IDs, ])
  rm(IDs)
  # check whether DFs have been synchronised
  stopifnot(all(sapply(2:4, function(x) identical(FILES[[1]][, 1:2], FILES[[x]][, 1:2]))))
  # write output
  for (i in 1:4) {
    write.table(FILES[[i]], file=names(FILES)[i], sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  }; rm(i)
}

#' Extract both logR and BAF values from sequencing data
#'
#' Method derived from the Battenberg package (https://github.com/Wedge-lab/battenberg).
#'
#' @param tumourseqfile Full path to the tumour BAM/CRAM file.
#' @param normalseqfile Full path to the normal BAM/CRAM file.
#' @param tumourname Identifier to be used for tumour output files.
#' @param normalname Identifier to be used for normal output files.
#' @param allelecounter_exe Path to the allele counter executable.
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr").
#' @param loci.prefix Prefix path to the loci data (e.g. "G1000_loci_chr").
#' @param gender Gender information, either 'XX' (=female) or 'XY' (=male).
#' @param genomeVersion Genome version, either 'hg19' or 'hg38'.
#' @param nthreads The number of parallel processes for getting allele counts (optional, default=1).
#' @param tumourLogR_file Path to the tumour logR output (optional, paste0(tumourname, "_tumourLogR.txt")).
#' @param tumourBAF_file Path to the tumour BAF output (optional, paste0(tumourname, "_tumourBAF.txt")).
#' @param normalLogR_file Path to the normal logR output (optional, paste0(tumourname, "_normalLogR.txt")).
#' @param normalBAF_file Path to the normal BAF output (optional, paste0(tumourname, "_normalBAF.txt")).
#' @param minCounts Minimum depth required in the normal for a SNP to be considered (optional, default=10).
#' @param BED_file A BED file for only looking at SNPs within specific intervals (optional, default=NA).
#' @param probloci_file A file (chromosome <tab> position; no header) containing specific loci to ignore (optional, default=NA).
#' @param chrom_names A vector containing the names of chromosomes to be considered (optional, default=c(1:22, 'X')).
#' @param min_base_qual Minimum base quality required for a read to be counted (optional, default=20).
#' @param min_map_qual Minimum mapping quality required for a read to be counted (optional, default=35).
#' @param ref.fasta FASTA file used for generating CRAMs (optional, default=NA).
#' @param skip_allele_counting_tumour Flag, set to TRUE if tumour allele counting is already complete (files are expected in the working directory on disk; optional, default=FALSE).
#' @param skip_allele_counting_normal Flag, set to TRUE if normal allele counting is already complete (files are expected in the working directory on disk; optional, default=FALSE).
#' @param seed A seed to be set when randomising the alleles (optional, default=as.integer(Sys.time())).
#' @author sd11, tl
#' @export
ascat.prepareHTS = function(tumourseqfile, normalseqfile, tumourname, normalname, allelecounter_exe, alleles.prefix, loci.prefix, gender, genomeVersion,
                            nthreads=1, tumourLogR_file=NA, tumourBAF_file=NA, normalLogR_file=NA, normalBAF_file=NA, minCounts=10, BED_file=NA,
                            probloci_file=NA, chrom_names=c(1:22, "X"), min_base_qual=20, min_map_qual=35, ref.fasta=NA,
                            skip_allele_counting_tumour=FALSE, skip_allele_counting_normal=FALSE, seed=as.integer(Sys.time())) {
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  registerDoParallel(cores=nthreads)

  if (is.na(tumourLogR_file)) tumourLogR_file=paste0(tumourname, "_tumourLogR.txt")
  if (is.na(tumourBAF_file)) tumourBAF_file=paste0(tumourname, "_tumourBAF.txt")
  if (is.na(normalLogR_file)) normalLogR_file=paste0(tumourname, "_normalLogR.txt")
  if (is.na(normalBAF_file)) normalBAF_file=paste0(tumourname, "_normalBAF.txt")

  if (!skip_allele_counting_tumour) {
    # Obtain allele counts at specific loci for tumour
    foreach(CHR=chrom_names) %dopar% {
      CHR=get("CHR")
      ascat.getAlleleCounts(seq.file=tumourseqfile,
                            output.file=paste0(tumourname, "_alleleFrequencies_chr", CHR, ".txt"),
                            loci.file=paste0(loci.prefix, CHR, ".txt"),
                            min.base.qual=min_base_qual,
                            min.map.qual=min_map_qual,
                            allelecounter.exe=allelecounter_exe,
                            ref.fasta=ref.fasta)
    }
  }
  if (!skip_allele_counting_normal) {
    # Obtain allele counts at specific loci for normal
    foreach(CHR=chrom_names) %dopar% {
      CHR=get("CHR")
      ascat.getAlleleCounts(seq.file=normalseqfile,
                            output.file=paste0(normalname, "_alleleFrequencies_chr", CHR, ".txt"),
                            loci.file=paste0(loci.prefix, CHR, ".txt"),
                            min.base.qual=min_base_qual,
                            min.map.qual=min_map_qual,
                            allelecounter.exe=allelecounter_exe,
                            ref.fasta=ref.fasta)
    }
  }
  # Obtain BAF and LogR from the raw allele counts
  ascat.getBAFsAndLogRs(samplename=tumourname,
                        tumourAlleleCountsFile.prefix=paste0(tumourname, "_alleleFrequencies_chr"),
                        normalAlleleCountsFile.prefix=paste0(normalname, "_alleleFrequencies_chr"),
                        tumourLogR_file=tumourLogR_file,
                        tumourBAF_file=tumourBAF_file,
                        normalLogR_file=normalLogR_file,
                        normalBAF_file=normalBAF_file,
                        alleles.prefix=alleles.prefix,
                        gender=gender,
                        genomeVersion=genomeVersion,
                        chrom_names=chrom_names,
                        minCounts=minCounts,
                        BED_file=BED_file,
                        probloci_file=probloci_file,
                        seed=seed)

  # Synchronise all information
  ascat.synchroniseFiles(samplename=tumourname,
                         tumourLogR_file=tumourLogR_file,
                         tumourBAF_file=tumourBAF_file,
                         normalLogR_file=normalLogR_file,
                         normalBAF_file=normalBAF_file)
}

#' Function to concatenate allele counter output
#' @noRd
readAlleleCountFiles=function(prefix, suffix, chrom_names, minCounts, keep_chr_string=FALSE) {
  files=paste0(prefix, chrom_names, suffix)
  files=files[sapply(files, function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind, lapply(files, function(x) {
    tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=TRUE), stringsAsFactors=FALSE)
    tmp=tmp[tmp[, 7]>=minCounts, ]
    if (nrow(tmp)>0) {
      if (!keep_chr_string) tmp[, 1]=gsub("^chr", "", tmp[, 1])
      rownames(tmp)=paste0(tmp[, 1], "_", tmp[, 2])
    }
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}

#' Function to concatenate all alleles files
#' @noRd
readAllelesFiles=function(prefix, suffix, chrom_names, add_chr_string=FALSE) {
  files=paste0(prefix, chrom_names, suffix)
  files=files[sapply(files, function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind, lapply(files, function(x) {
    tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=TRUE))
    tmp=tmp[!is.na(tmp[, 2] & !is.na(tmp[, 3])), ]
    tmp=tmp[!duplicated(tmp[, 1]), ]
    tmp$chromosome=gsub(paste0(prefix, "(", paste(chrom_names, collapse="|"), ")", suffix), "\\1", x)
    if (add_chr_string) tmp$chromosome=paste0("chr", tmp$chromosome)
    tmp=tmp[, c(4, 1:3)]
    rownames(tmp)=paste0(tmp[, 1], "_", tmp[, 2])
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}

#' Function to concatenate all loci files
#' @noRd
readLociFiles=function(prefix, suffix, chrom_names, keep_chr_string=FALSE) {
  files=paste0(prefix, chrom_names, suffix)
  files=files[sapply(files, function(x) file.exists(x) && file.info(x)$size>0)]
  stopifnot(length(files)>0)
  data=do.call(rbind, lapply(files, function(x) {
    tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=FALSE))
    if (!keep_chr_string) tmp[, 1]=gsub("^chr", "", tmp[, 1])
    rownames(tmp)=paste0(tmp[, 1], "_", tmp[, 2])
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}
