#' @title ascat.runAscat
#' @description ASCAT main function, calculating the allele-specific copy numbers
#' @param ASCATobj an ASCAT object from ascat.aspcf
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100\% aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
#' @param pdfPlot Optional flag if nonrounded plots and ASCAT profile in pdf format are desired. Default=F
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
# @param textFlag Optional flag to add the positions of fragments located outside of the plotting area to the plots. Default=F
#' @param circos Optional file to output the non-rounded values in Circos track format. Default=NA
#' @param min_ploidy optional numerical parameter determining the minimum boundary of the ploidy solution search space (expert parameter, don't adapt unless you know what you're doing). Default=1.5
#' @param max_ploidy optional numerical parameter determining the maximum boundary of the ploidy solution search space (expert parameter, don't adapt unless you know what you're doing). Default=5.5
#' @param min_purity optional numerical parameter determining the minimum boundary of the purity solution search space (expert parameter, don't adapt unless you know what you're doing). Default=0.1
#' @param max_purity optional numerical parameter determining the maximum boundary of the purity solution search space (expert parameter, don't adapt unless you know what you're doing). Default=1.05
#' @param rho_manual optional argument to override ASCAT optimization and supply rho parameter (expert parameter, don't adapt unless you know what you're doing).
#' @param psi_manual optional argument to override ASCAT optimization and supply psi parameter (expert parameter, don't adapt unless you know what you're doing).
#' @param img.dir directory in which figures will be written
#' @param img.prefix prefix for figure names
#' @param write_segments Optional flag to output segments in text files (.segments_raw.txt and .segments.txt under img.dir). Default=F
#' @details Note: for copy number only probes, nA contains the copy number value and nB = 0.
#' @return an ASCAT output object, containing:\cr
#' 1. nA: copy number of the A allele\cr
#' 2. nB: copy number of the B allele\cr
#' 3. purity: the tumour purity of all arrays\cr
#' 4. aberrantcellfraction: the aberrant cell fraction (=tumour purity) of all arrays\cr
#' 5. ploidy: the ploidy of all arrays\cr
#' 6. failedarrays: arrays on which ASCAT analysis failed\cr
#' 7. nonaberrantarrays: arrays on which ASCAT analysis indicates that they show virtually no aberrations\cr
#' 8. segments: an array containing the copy number segments of each sample (not including failed arrays)\cr
#' 9. segments_raw: an array containing the copy number segments of each sample without any rounding applied\cr
#' 10. distance_matrix: distances for a range of ploidy and tumor percentage values
#'
#' @export
#'
ascat.runAscat = function(ASCATobj, gamma = 0.55, pdfPlot = FALSE, y_limit = 5, circos=NA, min_ploidy=1.5, max_ploidy=5.5, min_purity=0.1, max_purity=1.05, rho_manual = NA, psi_manual = NA, img.dir=".", img.prefix="", write_segments=FALSE) {
  goodarrays=NULL
  N_samples=dim(ASCATobj$Tumor_LogR)[2]
  res = vector("list", N_samples)
  stopifnot(length(rho_manual)==length(psi_manual)) # check consistency
  if (length(rho_manual)==1 && is.na(rho_manual) && N_samples>1) {
    rho_manual=rep(NA, N_samples)
    psi_manual=rep(NA, N_samples)
  } else {
    stopifnot(length(rho_manual)==N_samples)
  }
  for (arraynr in 1:N_samples) {
    print.noquote(paste("Sample ", ASCATobj$samples[arraynr], " (", arraynr, "/", length(ASCATobj$samples), ")", sep=""))
    lrr=ASCATobj$Tumor_LogR[, arraynr]
    names(lrr)=rownames(ASCATobj$Tumor_LogR)
    baf=ASCATobj$Tumor_BAF[, arraynr]
    names(baf)=rownames(ASCATobj$Tumor_BAF)
    lrrsegm = ASCATobj$Tumor_LogR_segmented[, arraynr]
    names(lrrsegm) = rownames(ASCATobj$Tumor_LogR_segmented)
    bafsegm = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop=FALSE]
    names(bafsegm) = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    failedqualitycheck = FALSE
    if (ASCATobj$samples[arraynr] %in% ASCATobj$failedarrays) {
      failedqualitycheck = TRUE
    }
    ending = ifelse(pdfPlot, "pdf", "png")
    circosName=NA
    if (!is.na(circos)) {
      circosName=paste(circos, "_", ASCATobj$samples[arraynr], sep="")
    }
    res[[arraynr]] = runASCAT(lrr, baf, lrrsegm, bafsegm, ASCATobj$gender[arraynr], ASCATobj$SNPpos, ASCATobj$ch, ASCATobj$chrs, ASCATobj$sexchromosomes, failedqualitycheck,
                              file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".sunrise.png", sep="")), file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".ASCATprofile.", ending, sep="")),
                              file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".rawprofile.", ending, sep="")), NA,
                              gamma, rho_manual[arraynr], psi_manual[arraynr], pdfPlot, y_limit, circosName, min_ploidy, max_ploidy, min_purity, max_purity, ASCATobj$X_nonPAR)
    if (!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays)+1] = arraynr
    }
  }

  if (length(goodarrays)>0) {
    n1 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    n2 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) = rownames(ASCATobj$Tumor_LogR)
    rownames(n2) = rownames(ASCATobj$Tumor_LogR)
    colnames(n1) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    colnames(n2) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[, i] = res[[goodarrays[i]]]$nA
      n2[, i] = res[[goodarrays[i]]]$nB
    }

    distance_matrix = vector("list", length(goodarrays))
    names(distance_matrix) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      distance_matrix[[i]] = res[[goodarrays[i]]]$distance_matrix
    }

    tp = vector(length=length(goodarrays))
    psi = vector(length=length(goodarrays))
    ploidy = vector(length=length(goodarrays))
    goodnessOfFit = vector(length=length(goodarrays))
    naarrays = NULL
    for (i in 1:length(goodarrays)) {
      tp[i] = res[[goodarrays[i]]]$rho
      psi[i] = res[[goodarrays[i]]]$psi
      ploidy[i] = mean(res[[goodarrays[i]]]$nA+res[[goodarrays[i]]]$nB, na.rm=TRUE)
      goodnessOfFit[i] = res[[goodarrays[i]]]$goodnessOfFit
      if (res[[goodarrays[i]]]$nonaberrant) {
        naarrays = c(naarrays, ASCATobj$samples[goodarrays[i]])
      }
    }
    fa = colnames(ASCATobj$Tumor_LogR)[-goodarrays]
    names(tp) = colnames(n1)
    names(ploidy) = colnames(n1)
    names(psi) = colnames(n1)
    names(goodnessOfFit) = colnames(n1)

    seg = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg
      seg = rbind(seg, cbind(ASCATobj$samples[goodarrays[i]], as.vector(ASCATobj$SNPpos[segje[, 1], 1]),
                             ASCATobj$SNPpos[segje[, 1], 2],
                             ASCATobj$SNPpos[segje[, 2], 2], segje[, 3], segje[, 4]))
    }
    colnames(seg) = c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor")
    seg = data.frame(seg, stringsAsFactors=FALSE)
    seg[, 3]=as.numeric(seg[, 3])
    seg[, 4]=as.numeric(seg[, 4])
    seg[, 5]=as.numeric(seg[, 5])
    seg[, 6]=as.numeric(seg[, 6])

    seg_raw = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg_raw
      seg_raw = rbind(seg_raw, cbind(ASCATobj$samples[goodarrays[i]], as.vector(ASCATobj$SNPpos[segje[, 1], 1]),
                                     ASCATobj$SNPpos[segje[, 1], 2],
                                     ASCATobj$SNPpos[segje[, 2], 2], segje[, 3], segje[, 4:ncol(segje)]))

    }
    colnames(seg_raw) = c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor", "nAraw", "nBraw")

    seg_raw = data.frame(seg_raw, stringsAsFactors=FALSE)
    seg_raw[, 3]=as.numeric(seg_raw[, 3])
    seg_raw[, 4]=as.numeric(seg_raw[, 4])
    seg_raw[, 5]=as.numeric(seg_raw[, 5])
    seg_raw[, 6]=as.numeric(seg_raw[, 6])
    seg_raw[, 7]=as.numeric(seg_raw[, 7])
    seg_raw[, 8]=as.numeric(seg_raw[, 8])

    if (write_segments) {
      for (i in 1:length(goodarrays)) {
        # Write rounded segments
        mySeg=data.frame(sample=rep(ASCATobj$samples[goodarrays[i]], nrow(res[[goodarrays[i]]]$seg)),
                         chr=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg[, 1], 1],
                         startpos=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg[, 1], 2],
                         endpos=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg[, 2], 2],
                         nMajor=res[[goodarrays[i]]]$seg[, 3],
                         nMinor=res[[goodarrays[i]]]$seg[, 4],
                         stringsAsFactors=FALSE)
        write.table(mySeg, file=file.path(img.dir, paste0(ASCATobj$samples[goodarrays[i]], ".segments.txt")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        # Write unrounded segments
        mySeg=data.frame(sample=rep(ASCATobj$samples[goodarrays[i]], nrow(res[[goodarrays[i]]]$seg_raw)),
                         chr=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg_raw[, 1], 1],
                         startpos=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg_raw[, 1], 2],
                         endpos=ASCATobj$SNPpos[res[[goodarrays[i]]]$seg_raw[, 2], 2],
                         nMajor=res[[goodarrays[i]]]$seg_raw[, 3],
                         nMinor=res[[goodarrays[i]]]$seg_raw[, 4],
                         nAraw=res[[goodarrays[i]]]$seg_raw[, 5],
                         nBraw=res[[goodarrays[i]]]$seg_raw[, 6],
                         stringsAsFactors=FALSE)
        write.table(mySeg, file=file.path(img.dir, paste0(ASCATobj$samples[goodarrays[i]], ".segments_raw.txt")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        rm(mySeg)
      }; rm(i)
    }

  } else {
    n1 = NULL
    n2 = NULL
    tp = NULL
    ploidy = NULL
    psi = NULL
    goodnessOfFit = NULL
    fa = colnames(ASCATobj$Tumor_LogR)
    naarrays = NULL
    seg = NULL
    seg_raw = NULL
    distance_matrix = NULL
  }

  return(list(nA = n1, nB = n2, purity = tp, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, goodnessOfFit = goodnessOfFit,
              failedarrays = fa, nonaberrantarrays = naarrays, segments = seg, segments_raw = seg_raw, distance_matrix = distance_matrix))
}

#' @title runASCAT
#' @description the ASCAT main function
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param gender a vector of gender for each cases ("XX" or "XY"). Default = NULL: all female ("XX")
#' @param SNPpos position of all SNPs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param chrnames a vector containing the names for the chromosomes (e.g. c(1:22, "X"))
#' @param sexchromosomes a vector containing the names for the sex chromosomes
#' @param failedqualitycheck did the sample fail any previous quality check or not?
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file
#' @param aberrationreliabilitypng aberration reliability score is plotted if filename is given
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100\% aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
#' @param rho_manual optional argument to override ASCAT optimization and supply rho parameter (not recommended)
#' @param psi_manual optional argument to override ASCAT optimization and supply psi parameter (not recommended)
#' @param pdfPlot Optional flag if nonrounded plots and ASCAT profile in pdf format are desired. Default=FALSE
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#' @param circos Optional file to output the non-rounded values in Circos track format. Default=NA
#' @param min_ploidy a numerical parameter determining the minimum boundary of the ploidy solution search space. Default=1.5
#' @param max_ploidy a numerical parameter determining the maximum boundary of the ploidy solution search space. Default=5.5
#' @param min_purity a numerical parameter determining the minimum boundary of the purity solution search space. Default=0.1
#' @param max_purity a numerical parameter determining the maximum boundary of the purity solution search space. Default=1.05
#' @param X_nonPAR Optional vector containing genomic coordinates (start & stop) of nonPAR region on X. Default=NULL
#'
#' @keywords internal
#'
#' @return list containing optimal purity and ploidy
#'
#' @import RColorBrewer
#'
#' @export
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames, sexchromosomes, failedqualitycheck = FALSE,
                    distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, aberrationreliabilitypng = NA, gamma = 0.55,
                    rho_manual = NA, psi_manual = NA, pdfPlot = FALSE, y_limit = 5, circos=NA, min_ploidy=1.5, max_ploidy=5.5, min_purity=0.1,
                    max_purity=1.05, X_nonPAR=NULL) {
  ch = chromosomes
  chrs = chrnames
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  SNPposhet = SNPpos[names(bafsegmented), ]
  autoprobes = !(SNPposhet[, 1] %in% sexchromosomes)

  b2 = b[autoprobes]
  r2 = r[autoprobes]

  s = make_segments(r2, b2)
  d = create_distance_matrix(s, gamma, min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity)
  plot_d=d

  TheoretMaxdist = sum(rep(0.25, dim(s)[1]) * s[, "length"] * ifelse(s[, "b"]==0.5, 0.05, 1), na.rm=TRUE)

  # flag the sample as non-aberrant if necessary
  nonaberrant = FALSE
  MINABB = 0.03
  MINABBREGION = 0.005

  percentAbb = sum(ifelse(s[, "b"]==0.5, 0, 1)*s[, "length"])/sum(s[, "length"])
  maxsegAbb = max(ifelse(s[, "b"]==0.5, 0, s[, "length"]))/sum(s[, "length"])
  if (percentAbb <= MINABB && maxsegAbb <= MINABBREGION) {
    nonaberrant = TRUE
  }


  MINPLOIDY = min_ploidy
  MAXPLOIDY = max_ploidy
  MINRHO = 0.2
  MINGOODNESSOFFIT = 80
  MINPERCZERO = 0.02
  MINPERCZEROABB = 0.1
  MINPERCODDEVEN = 0.05
  MINPLOIDYSTRICT = 1.7
  MAXPLOIDYSTRICT = 2.3

  nropt = 0
  localmin = NULL
  optima = list()

  if (!failedqualitycheck && is.na(rho_manual)) {

    # first, try with all filters
    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i, j]
        seld = d[(i-3):(i+3), (j-3):(j+3)]
        seld[4, 4] = max(seld)
        if (min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])
          nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
          nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho

          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])

          percentzero = (sum((round(nA)==0)*s[, "length"])+sum((round(nB)==0)*s[, "length"]))/sum(s[, "length"])

          goodnessOfFit = (1-m/TheoretMaxdist) * 100

          if (!nonaberrant && ploidy > MINPLOIDY && ploidy < MAXPLOIDY && rho >= MINRHO && goodnessOfFit > MINGOODNESSOFFIT && percentzero > MINPERCZERO) {
            nropt = nropt + 1
            optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }

    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if (nropt == 0  && MINPLOIDY < MAXPLOIDYSTRICT && MAXPLOIDY > MINPLOIDYSTRICT) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i, j]
          seld = d[(i-3):(i+3), (j-3):(j+3)]
          seld[4, 4] = max(seld)
          if (min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
            nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])

            perczeroAbb = (sum((round(nA)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))+sum((round(nB)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1)))/sum(s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT && ploidy < MAXPLOIDYSTRICT && rho >= MINRHO && goodnessOfFit > MINGOODNESSOFFIT && perczeroAbb > MINPERCZEROABB) {
              nropt = nropt + 1
              optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }

    # if still no solution, allow solutions with 100% aberrant cells (include the borders with rho = 1), but in first instance, keep the percentzero > 0.01 filter
    if (nropt == 0) {
      #first, include borders
      cold = which(as.numeric(colnames(d))>1)
      d[, cold]=1E20
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i, j]
          seld = d[(i-3):(i+3), (j-3):(j+3)]
          seld[4, 4] = max(seld)
          if (min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
            nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])

            percentzero = (sum((round(nA)==0)*s[, "length"])+sum((round(nB)==0)*s[, "length"]))/sum(s[, "length"])
            percOddEven = sum((round(nA) %% 2 == 0 & round(nB) %% 2 == 1 | round(nA) %% 2 == 1 & round(nB) %% 2 == 0)*s[, "length"])/sum(s[, "length"])
            perczeroAbb = (sum((round(nA)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))+sum((round(nB)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1)))/sum(s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (!nonaberrant && ploidy > MINPLOIDY && ploidy < MAXPLOIDY && rho >= MINRHO && goodnessOfFit > MINGOODNESSOFFIT &&
                  (perczeroAbb > MINPERCZEROABB || percentzero > MINPERCZERO || percOddEven > MINPERCODDEVEN)) {
              nropt = nropt + 1
              optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }

    # if still no solution, drop the percentzero > MINPERCENTZERO filter, but strict ploidy borders
    if (nropt == 0  && MINPLOIDY < MAXPLOIDYSTRICT && MAXPLOIDY > MINPLOIDYSTRICT) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i, j]
          seld = d[(i-3):(i+3), (j-3):(j+3)]
          seld[4, 4] = max(seld)
          if (min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
            nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])

            perczeroAbb = (sum((round(nA)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))+sum((round(nB)==0)*s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1)))/sum(s[, "length"]*ifelse(s[, "b"]==0.5, 0, 1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT && ploidy < MAXPLOIDYSTRICT && rho >= MINRHO && goodnessOfFit > MINGOODNESSOFFIT) {
              nropt = nropt + 1
              optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }
  }

  if (!is.na(rho_manual)) {

    rho = rho_manual
    psi = psi_manual

    nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
    nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho

    # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
    ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])

    nMinor = NULL
    if (sum(nA, na.rm=TRUE) < sum(nB, na.rm=TRUE)) {
      nMinor = nA
    } else {
      nMinor = nB
    }
    m = sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * ifelse(s[, "b"]==0.5, 0.05, 1), na.rm=TRUE)

    goodnessOfFit = (1-m/TheoretMaxdist) * 100

    nropt = 1
    optima[[1]] = c(m, rho, psi, ploidy, goodnessOfFit)
    localmin[1] = m

  }


  if (nropt>0) {
    if (is.na(rho_manual)) {
      optlim = sort(localmin)[1]
      for (i in 1:length(optima)) {
        if (optima[[i]][1] == optlim) {
          psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
          rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
          if (rho_opt1 > 1) {
            rho_opt1 = 1
          }
          #ploidy_opt1 = optima[[i]][4]
          goodnessOfFit_opt1 = optima[[i]][5]
        }
      }
    } else {
      rho_opt1 = optima[[1]][2]
      psi_opt1 = optima[[1]][3]
      #ploidy_opt1 = optima[[1]][4]
      goodnessOfFit_opt1 = optima[[1]][5]
    }
  }

  if (nropt>0) {
    #plot Sunrise
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    }
    ascat.plotSunrise(plot_d, psi_opt1, rho_opt1)
    if (!is.na(distancepng)) {
      dev.off()
    }

    rho = rho_opt1
    psi = psi_opt1
    SNPposhet = SNPpos[names(bafsegmented), ]
    haploidchrs = unique(c(substring(gender, 1, 1), substring(gender, 2, 2)))
    if (substring(gender, 1, 1)==substring(gender, 2, 2)) {
      haploidchrs = setdiff(haploidchrs, substring(gender, 1, 1))
    }
    diploidprobes = !(SNPposhet[, 1] %in% haploidchrs)
    if (!is.null(X_nonPAR) && gender=="XY") diploidprobes=diploidprobes_fixnonPAR(diploidprobes, SNPposhet, X_nonPAR, paste0(r, "/", b[, 1]))
    nullchrs = setdiff(sexchromosomes, unique(c(substring(gender, 1, 1), substring(gender, 2, 2))))
    nullprobes = SNPposhet[, 1] %in% nullchrs

    nAfull = ifelse(diploidprobes,
                    (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes, 0,
                           ifelse(b<0.5, (rho-1 + ((1-rho)*2+rho*psi)*2^(r/gamma))/rho, 0)))
    nBfull = ifelse(diploidprobes,
                    (rho-1+b*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes, 0,
                           ifelse(b<0.5, 0, (rho-1 + ((1-rho)*2+rho*psi)*2^(r/gamma))/rho)))
    nA = pmax(round(nAfull), 0)
    nB = pmax(round(nBfull), 0)

    if (!is.na(circos)) {
      frame<-cbind(SNPposhet, nAfull, nBfull)
      chrSegmA<-rle(frame$nAfull)
      chrSegmB<-rle(frame$nBfull)
      if (all(chrSegmA$lengths==chrSegmB$lengths)) {
        start=1
        for (i in 1:length(chrSegmA$values)) {
          valA<-chrSegmA$values[i]
          valB<-chrSegmB$values[i]
          size<-chrSegmA$lengths[i]
          write(c(paste("hs", frame[start, 1], sep=""), frame[start, 2], frame[(start+size-1), 2], valA), file = paste(circos, "_major", sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          write(c(paste("hs", frame[start, 1], sep=""), frame[start, 2], frame[(start+size-1), 2], valB), file = paste(circos, "_minor", sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          start=start+size
        }
      } else {
        print("Major and minor allele copy numbers are segmented differently.")
      }
    }

    rho = rho_opt1
    psi = psi_opt1

    diploidprobes = !(SNPpos[, 1] %in% haploidchrs)
    if (!is.null(X_nonPAR) && gender=="XY") diploidprobes=diploidprobes_fixnonPAR(diploidprobes, SNPpos, X_nonPAR, lrrsegmented)
    nullprobes = SNPpos[, 1] %in% nullchrs

    #this replaces an occurrence of unique that caused problems
    #introduces segment spanning over chr ends, when two consecutive probes from diff chr have same logR!
    # build helping vector
    chrhelp = vector(length=length(lrrsegmented))
    for (chrnr in 1:length(ch)) {
      chrke = ch[[chrnr]]
      chrhelp[chrke] = chrnr
    }

    tlr2 = rle(lrrsegmented)
    tlr.chr= rle(chrhelp)

    tlrstart = c(1, cumsum(tlr2$lengths)+1)
    tlrstart = tlrstart[1:(length(tlrstart)-1)]
    tlrend = cumsum(tlr2$lengths)

    tlrstart.chr= c(1, cumsum(tlr.chr$lengths)+1)
    tlrstart.chr = tlrstart.chr[1:(length(tlrstart.chr)-1)]
    tlrend.chr = cumsum(tlr.chr$lengths)

    tlrend<-sort(union(tlrend, tlrend.chr))
    tlrstart<-sort(union(tlrstart, tlrstart.chr))

    tlr=NULL
    for (ind in tlrstart) {
      val<-lrrsegmented[ind]
      tlr<-c(tlr, val)
    }

    # For each LRR probe, find the matching BAF probe
    # and its position in bafsegmented
    probeLookup = data.frame(
      lrrprobe = names(lrrsegmented),
      bafpos = match(names(lrrsegmented), names(bafsegmented)),
      stringsAsFactors=FALSE
    )

    seg = NULL
    for (i in 1:length(tlr)) {
      logR = tlr[i]
      #pr = which(lrrsegmented==logR) # this was a problem
      pr = tlrstart[i]:tlrend[i]
      start = min(pr)
      end = max(pr)

      bafpos = probeLookup$bafpos[pr]
      bafpos = bafpos[!is.na(bafpos)]
      bafke  = bafsegmented[bafpos][1]

      #if bafke is NA, this means that we are dealing with a germline homozygous stretch with a copy number change within it.
      #in this case, nA and nB are irrelevant, just their sum matters
      if (is.na(bafke)) {
        bafke=0
      }

      nAraw = ifelse(diploidprobes[start],
                     (rho-1 - (bafke-1)*2^(logR/gamma) * ((1-rho)*2+rho*psi)) / rho,
                     ifelse(nullprobes[start], 0,
                            (rho-1 + ((1-rho)*2+rho*psi)*2^(logR/gamma))/rho))
      nBraw = ifelse(diploidprobes[start], (rho-1+bafke*2^(logR/gamma) * ((1-rho)*2+rho*psi))/rho, 0)
      # correct for negative values:
      if (nAraw+nBraw<0) {
        nAraw = 0
        nBraw = 0
      } else if (nAraw<0) {
        nBraw = nAraw+nBraw
        nAraw = 0
      } else if (nBraw<0) {
        nAraw = nAraw+nBraw
        nBraw = 0
      }
      # when evidence for odd copy number in segments of BAF = 0.5, assume a deviation..
      limitround = 0.5
      nA = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nAraw)+1,
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nAraw),
                                round(nAraw))),
                  round(nAraw))
      nB = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nBraw),
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nBraw)-1,
                                round(nBraw))),
                  round(nBraw))
      if (is.null(seg)) {
        seg = t(as.matrix(c(start, end, nA, nB)))
        seg_raw = t(as.matrix(c(start, end, nA, nB, nAraw, nBraw)))
      } else {
        seg = rbind(seg, c(start, end, nA, nB))
        seg_raw = rbind(seg_raw, c(start, end, nA, nB, nAraw, nBraw))
      }
    }
    colnames(seg)=c("start", "end", "nA", "nB")
    colnames(seg_raw)=c("start", "end", "nA", "nB", "nAraw", "nBraw")

    # every repeat joins 2 ends. 20 repeats will join about 1 million ends..
    for (rep in 1:20) {
      seg2=seg
      seg = NULL
      skipnext = FALSE
      for (i in 1:dim(seg2)[1]) {
        if (!skipnext) {
          if (i != dim(seg2)[1] && seg2[i, "nA"]==seg2[i+1, "nA"] && seg2[i, "nB"]==seg2[i+1, "nB"] &&
                chrhelp[seg2[i, "end"]]==chrhelp[seg2[i+1, "start"]]) {
            segline = c(seg2[i, "start"], seg2[i+1, "end"], seg2[i, 3:4])
            skipnext = TRUE
          } else {
            segline = seg2[i, ]
          }

          if (is.null(seg)) {
            seg = t(as.matrix(segline))
          } else {
            seg = rbind(seg, segline)
          }
        } else {
          skipnext = FALSE
        }
      }
      colnames(seg)=colnames(seg2)
    }
    rownames(seg)=NULL

    nMajor = vector(length = length(lrrsegmented))
    names(nMajor) = names(lrrsegmented)
    nMinor = vector(length = length(lrrsegmented))
    names(nMinor) = names(lrrsegmented)

    for (i in 1:dim(seg)[1]) {
      nMajor[seg[i, "start"]:seg[i, "end"]] = seg[i, "nA"]
      nMinor[seg[i, "start"]:seg[i, "end"]] = seg[i, "nB"]
    }

    n1all = vector(length = length(lrrsegmented))
    names(n1all) = names(lrrsegmented)
    n2all = vector(length = length(lrrsegmented))
    names(n2all) = names(lrrsegmented)

    # note: any of these can have length 0
    NAprobes = which(is.na(lrr))
    heteroprobes = setdiff(which(names(lrrsegmented) %in% names(bafsegmented)), NAprobes)
    homoprobes = setdiff(setdiff(which(!is.na(baf)), heteroprobes), NAprobes)
    CNprobes = setdiff(which(is.na(baf)), NAprobes)

    n1all[NAprobes] = NA
    n2all[NAprobes] = NA
    n1all[CNprobes] = nMajor[CNprobes]+nMinor[CNprobes]
    n2all[CNprobes] = 0
    heteroprobes2 = names(lrrsegmented)[heteroprobes]
    n1all[heteroprobes] = ifelse(baf[heteroprobes2]<=0.5, nMajor[heteroprobes], nMinor[heteroprobes])
    n2all[heteroprobes] = ifelse(baf[heteroprobes2]>0.5, nMajor[heteroprobes], nMinor[heteroprobes])
    n1all[homoprobes] = ifelse(baf[homoprobes]<=0.5, nMajor[homoprobes]+nMinor[homoprobes], 0)
    n2all[homoprobes] = ifelse(baf[homoprobes]>0.5, nMajor[homoprobes]+nMinor[homoprobes], 0)

    # plot nonrounded ASCAT profile
    if (is.na(nonroundedprofilepng)) {
      dev.new(10, 5)
    } else {
      if (pdfPlot) {
        pdf(file = nonroundedprofilepng, width = 20, height = y_limit, pointsize=20)
      } else {
        png(filename = nonroundedprofilepng, width = 2000, height = (y_limit*100), res = 200)
      }
    }
    ascat.plotNonRounded(mean(n1all+n2all, na.rm=TRUE), rho_opt1, goodnessOfFit_opt1, nonaberrant, nAfull, nBfull, y_limit, bafsegmented, ch, lrr, chrnames)

    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }

    # plot ASCAT profile
    if (is.na(copynumberprofilespng)) {
      dev.new(10, 2.5)
    } else {
      if (pdfPlot) {
        pdf(file = copynumberprofilespng, width = 20, height = y_limit, pointsize=20)
      } else {
        png(filename = copynumberprofilespng, width = 2000, height = (y_limit*100), res = 200)
      }
    }
    #plot ascat profile
    ascat.plotAscatProfile(n1all, n2all, heteroprobes, mean(n1all+n2all, na.rm=TRUE), rho_opt1, goodnessOfFit_opt1, nonaberrant, y_limit, ch, lrr, bafsegmented, chrnames)


    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }


    if (!is.na(aberrationreliabilitypng)) {
      png(filename = aberrationreliabilitypng, width = 2000, height = 500, res = 200)
      par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)

      diploidprobes = !(SNPposhet[, 1] %in% haploidchrs)
      nullprobes = SNPposhet[, 1] %in% nullchrs

      rBacktransform = ifelse(diploidprobes,
                              gamma*log((rho * (nA+nB) + (1-rho)*2) / ((1-rho)*2+rho*psi), 2),
                              # the value for nullprobes is arbitrary (but doesn't matter, as these are not plotted anyway because BAF=0.5)
                              ifelse(nullprobes, -10, gamma*log((rho * (nA+nB) + (1-rho)) / ((1-rho)*2+rho*psi), 2)))

      bBacktransform = ifelse(diploidprobes,
                              (1-rho+rho*nB) / (2-2*rho+ rho * (nA+nB)),
                              ifelse(nullprobes, 0.5, 0))

      rConf = ifelse(abs(rBacktransform)>0.15, pmin(100, pmax(0, 100 * (1-abs(rBacktransform-r)/abs(r)))), NA)
      bConf = ifelse(diploidprobes & bBacktransform!=0.5, pmin(100, pmax(0, ifelse(b==0.5, 100, 100 * (1-abs(bBacktransform-b)/abs(b-0.5))))), NA)
      confidence = ifelse(is.na(rConf), bConf, ifelse(is.na(bConf), rConf, (rConf+bConf)/2))
      maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f", mean(confidence, na.rm=TRUE)), "%", sep="")
      plot(c(1, length(nAfull)), c(0, 100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
      points(confidence, col="blue", pch = "|")
      abline(v=0, lty=1, col="lightgrey")
      chrk_tot_len = 0
      for (i in 1:length(ch)) {
        chrk = ch[[i]]
        chrk_hetero = intersect(names(lrr)[chrk], names(bafsegmented))
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk_hetero)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2
        text(tpos, 5, chrs[i], pos = 1, cex = 2)
        abline(v=vpos, lty=1, col="lightgrey")
      }
      dev.off()
    }

    return(list(rho = rho_opt1, psi = psi_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = nonaberrant,
                nA = n1all, nB = n2all, seg = seg, seg_raw = seg_raw, distance_matrix = d))

  } else {

    name=gsub(".sunrise.png", "", basename(distancepng))

    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    ascat.plotSunrise(plot_d, 0, 0)
    dev.off()

    warning(paste("ASCAT could not find an optimal ploidy and purity value for sample ", name, ".\n", sep=""))
    return(list(rho = NA, psi = NA, goodnessOfFit = NA, nonaberrant = FALSE, nA = NA, nB = NA, seg = NA, seg_raw = NA, distance_matrix = NA))
  }

}

#' @title make_segments
#' @description Function to make segments of constant LRR and BAF.\cr
#' This function is more general and does not depend on specific ASPCF output, it can also handle segmention performed on LRR and BAF separately
#' @param r segmented logR
#' @param b segmented BAF
#' @return segments of constant logR and BAF including their lengths
#' @keywords internal
#' @export
make_segments = function(r, b) {
  m = matrix(ncol = 2, nrow = length(b))
  m[, 1] = r
  m[, 2] = b
  m = as.matrix(na.omit(m))
  pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) = c("r", "b", "length")
  index = 0
  previousb = -1
  previousr = 1E10
  for (i in 1:dim(m)[1]) {
    if (m[i, 2] != previousb || m[i, 1] != previousr) {
      index=index+1
      count=1
      pcf_segments[index, "r"] = m[i, 1]
      pcf_segments[index, "b"] = m[i, 2]
    } else {
      count = count + 1
    }
    pcf_segments[index, "length"] = count
    previousb = m[i, 2]
    previousr = m[i, 1]
  }
  pcf_segments = as.matrix(na.omit(pcf_segments))[, , drop=FALSE]
  return(pcf_segments)
}

# function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
# input: segmented LRR and BAF and the value for gamma
create_distance_matrix = function(segments, gamma, min_ploidy=NULL, max_ploidy=NULL, min_purity=NULL, max_purity=NULL) {
  s = segments
  # get ploidy boundaries
  if (is.null(min_ploidy) || is.null(max_ploidy)) {
    psi_pos = seq(1, 6, 0.05)
  } else {
    psi_pos = seq(min_ploidy-0.5, max_ploidy+0.5, 0.05)
  }
  # get purity boundaries
  if (is.null(min_purity) || is.null(max_purity)) {
    rho_pos = seq(0.1, 1.05, 0.01)
  } else {
    rho_pos = seq(round(min_purity, 2), round(max_purity, 2), 0.01)
  }

  # set up distance matrix
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  #dmin = 1E20
  # get distances for each ploidy and purity combination
  for (i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for (j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      nA = (rho-1 - (s[, "b"]-1)*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
      nB = (rho-1 + s[, "b"]*2^(s[, "r"]/gamma) * ((1-rho)*2+rho*psi))/rho
      # choose the minor allele
      nMinor = NULL
      if (sum(nA, na.rm=TRUE) < sum(nB, na.rm=TRUE)) {
        nMinor = nA
      } else {
        nMinor = nB
      }
      d[i, j] = sum(abs(nMinor - pmax(round(nMinor), 0))^2 * s[, "length"] * ifelse(s[, "b"]==0.5, 0.05, 1), na.rm=TRUE)
    }
  }
  return(d)
}

#' Function to fix diploidprobes for X based on nonPAR and segmentation (males only).
#' Setting diploidprobes to TRUE for PAR regions would generate an issue whenever a segment spans both PAR and nonPAR regions.
#' This is because probes, within the same segment (so with the same logR/BAF information) will use different CN equations (driven by diploidprobes).
#' We use a conservative approach where we keep segments as they are and measure overlap with nonPAR.
#' If the overlap is big enough (>50% of the segment), the whole segment is assigned to nonPAR and diploidprobes=FALSE for all probes within that segment.
#' @noRd
diploidprobes_fixnonPAR=function(diploidprobes, SNP_data, nonPAR, seg_info) {
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  stopifnot(length(diploidprobes)==nrow(SNP_data) && length(diploidprobes)==length(seg_info))
  SNP_data[, 1]=gsub("^chr", "", SNP_data[, 1]) # remove "chr"
  if (!"X" %in% SNP_data[, 1]) return(diploidprobes)
  # First, set all SNPs to diploidprobes=TRUE for X
  diploidprobes[which(SNP_data[, 1]=="X")]=TRUE
  # Create a GRanges object with nonPAR information
  nonPAR=GRanges(seqnames="X", ranges=IRanges(start=nonPAR[1], end=nonPAR[2]))
  # Create a DF with chr (X only), pos and seg_info (logR/BAF or logR only)
  DATA=SNP_data
  DATA$segment=seg_info
  DATA=DATA[which(DATA[, 1]=="X"), ]
  # Run rle on seg_info to get segments
  RLE=cumsum(c(1, rle(DATA$segment)$lengths))
  # Create a DF where each row corresponds to one segment
  SEGMENTS=do.call(rbind, lapply(1:(length(RLE)-1), function(x) {
    return(data.frame(Start=DATA[, 2][RLE[x]],
                      End=DATA[, 2][RLE[x+1]-1],
                      Segment=DATA$segment[RLE[x]]))
  }))
  # Convert DF to GRanges
  SEGMENTS=GRanges(seqnames="X", ranges=IRanges(start=SEGMENTS$Start, end=SEGMENTS$End), Segment=SEGMENTS$Segment)
  # Identify segments mapping to nonPAR
  FO=findOverlaps(nonPAR, SEGMENTS)
  if (length(FO)==0) return(diploidprobes)
  # For each segment, compute the size of the overlap relative to segment size. Only keep segments with high overlap (>50%) with nonPAR
  FO=FO[width(pintersect(nonPAR[queryHits(FO)], SEGMENTS[subjectHits(FO)]))/width(SEGMENTS[subjectHits(FO)])>0.5]
  if (length(FO)==0) return(diploidprobes)
  SEGMENTS=data.frame(SEGMENTS[subjectHits(FO)])
  # For each such segment, set all SNPs within that segment to diploidprobes=FALSE
  for (i in 1:nrow(SEGMENTS)) {
    diploidprobes[which(SNP_data[, 1]=="X" & SNP_data[, 2]>=SEGMENTS$start[i] & SNP_data[, 2]<=SEGMENTS$end[i])]=FALSE
  }; rm(i)
  return(diploidprobes)
}
