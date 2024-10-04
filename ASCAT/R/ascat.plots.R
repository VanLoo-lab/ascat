#' @title ascat.plotRawData
#' @description Plots SNP array data
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' @param img.dir directory in which figures will be written
#' @param img.prefix prefix for figure names
#' @param logr.y_values define Y min and max values for logR track (optional; default: c(-2, 2))
#'
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#'
#' @export
ascat.plotRawData = function(ASCATobj, img.dir=".", img.prefix="", logr.y_values=c(-2, 2)) {
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[i], ".tumour.png", sep="")), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    plot(c(1, dim(ASCATobj$Tumor_LogR)[1]), logr.y_values, type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[, i], col="red")
    points(ASCATobj$Tumor_LogR[, i], col="#77000011")
    abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos, logr.y_values[2], ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5, lty=1, col="lightgrey")
    }
    plot(c(1, dim(ASCATobj$Tumor_BAF)[1]), c(0, 1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[, i], col="red")
    points(ASCATobj$Tumor_BAF[, i], col="#77000011")
    abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5, lty=1, col="lightgrey")
    }
    dev.off()
  }

  if (!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[i], ".germline.png", sep="")), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
      plot(c(1, dim(ASCATobj$Germline_LogR)[1]), logr.y_values, type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_LogR[, i], col="red")
      points(ASCATobj$Germline_LogR[, i], col="#77000011")
      abline(v=0.5, lty=1, col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]]
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2
        text(tpos, logr.y_values[2], ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5, lty=1, col="lightgrey")
      }
      plot(c(1, dim(ASCATobj$Germline_BAF)[1]), c(0, 1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_BAF[, i], col="red")
      points(ASCATobj$Germline_BAF[, i], col="#77000011")
      abline(v=0.5, lty=1, col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]]
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2
        text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5, lty=1, col="lightgrey")
      }
      dev.off()
    }
  }
}

#' @title ascat.plotSegmentedData
#' @description plots the SNP array data before and after segmentation
#'
#' @param ASCATobj an ASCAT object (e.g. from ascat.aspcf)
#' @param img.dir directory in which figures will be written
#' @param img.prefix prefix for figure names
#' @param logr.y_values define Y min and max values for logR track (optional; default: c(-2, 2))
#'
#' @return png files showing raw and segmented tumour logR and BAF
#'
#' @export
#'
ascat.plotSegmentedData = function(ASCATobj, img.dir=".", img.prefix="", logr.y_values=c(-2, 2)) {
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Select_nonNAs = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) = rownames(ASCATobj$Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".ASPCF.png", sep="")), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main=3, cex.axis = 2)
    r = ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr]
    beta = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop=FALSE]
    plot(c(1, length(r)), logr.y_values, type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", LogR", sep=""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "#77000011", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    points(r, col="#1b38ae", pch=16)
    abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos, logr.y_values[2], ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5, lty=1, col="lightgrey")
    }
    plot(c(1, length(beta)), c(0, 1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", BAF", sep=""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "#77000011", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    points(beta, col = "#1b38ae", pch=16)
    points(1-beta, col = "#1b38ae", pch=16)
    abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5, lty=1, col="lightgrey")
    }
    dev.off()
  }
}

#' ascat.plotSunrise
#'
#' @param d distance matrix for a range of ploidy and tumour percentage values
#' @param psi_opt1 optimal ploidy
#' @param rho_opt1 optimal purity
#' @param minim when set to true, optimal regions in the sunrise plot are depicted in blue; if set to false, colours are inverted and red corresponds to optimal values (default: TRUE)
#'
#' @return plot visualising range of ploidy and tumour percentage values
#' @export
ascat.plotSunrise<-function(d, psi_opt1, rho_opt1, minim=TRUE) {

  par(mar = c(5, 5, 0.5, 0.5), cex=0.75, cex.lab=2, cex.axis=2)

  if (minim) {
    hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  } else {
    hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
  }
  image(log(d), col = hmcol, axes = FALSE, xlab = "Ploidy", ylab = "Purity")

  ploidy_min<-as.numeric(rownames(d)[1])
  ploidy_max<-as.numeric(rownames(d)[nrow(d)])
  purity_min<-as.numeric(colnames(d)[1])
  purity_max<-as.numeric(colnames(d)[ncol(d)])

  PLOIDY_LABELS=seq(ploidy_min, ploidy_max, length.out=max(c((ploidy_max-ploidy_min)+1, 2)))
  axis(1, at = (PLOIDY_LABELS-ploidy_min) / (ploidy_max-ploidy_min), labels = round(PLOIDY_LABELS, 2))
  if (purity_min==0.1 && purity_max==1.05) {
    PURITY_LABELS=seq(purity_min, purity_max, by = 0.3) # default Y-axis
  } else {
    PURITY_LABELS=seq(purity_min, purity_max, length.out=4)
  }
  axis(2, at = (PURITY_LABELS-purity_min) / (purity_max-purity_min), labels = PURITY_LABELS)

  if (psi_opt1>0 && rho_opt1>0) {
    points((psi_opt1-ploidy_min) / (ploidy_max-ploidy_min), (rho_opt1-purity_min) / (purity_max-purity_min), col="green", pch=4, cex = 2)
  }
}

#' @title ascat.plotNonRounded
#'
#' @description Function plotting the unrounded ASCAT copy number over all chromosomes
#' @param ploidy ploidy of the sample
#' @param rho purity of the sample
#' @param goodnessOfFit estimated goodness of fit
#' @param nonaberrant boolean flag denoting non-aberrated samples
#' @param nAfull copy number major allele
#' @param nBfull copy number minor allele
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
# @param textFlag Optional flag to add the positions of fragments located outside of the plotting area to the plots. Default=FALSE
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22, "X"))
#' @param ch a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#'
#' @return plot showing the nonrounded copy number profile, using base plotting function
#'
#' @export
ascat.plotNonRounded <- function(ploidy, rho, goodnessOfFit, nonaberrant, nAfull, nBfull, y_limit=5, bafsegmented, ch, lrr, chrs) {
  maintitle = paste("Ploidy: ", sprintf("%1.2f", ploidy), ", purity: ", sprintf("%2.0f", rho*100), "%, goodness of fit: ", sprintf("%2.1f", goodnessOfFit), "%", ifelse(nonaberrant, ", non-aberrant", ""), sep="")
  nBfullPlot<-ifelse(nBfull<y_limit, nBfull, y_limit+0.1)
  nAfullPlot<-ifelse((nAfull+nBfull)<y_limit, nAfull+nBfull, y_limit+0.1)
  colourTotal = "#943CC3" # purple
  colourMinor = "#60AF36" # green
  base.gw.plot(bafsegmented, nAfullPlot, nBfullPlot, colourTotal, colourMinor, maintitle, ch, lrr, chrs, y_limit, twoColours=TRUE)
}

#  @title create.bb.plot.average
#
#  @param BAFvals B Allele Frequency for every probe with position
#  @param subclones Segments with chromosomal locations
#  @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#  @param ploidy ploidy of the sample
#  @param rho purity of the sample
#  @param goodnessOfFit estimated goodness of fit
#  @param pos_min
#  @param pos_max
#  @param segment_states_min Vector containing copy number per segment minor allele
#  @param segment_states_tot Vector containing copy number per segment total copy number
#  @param chr.segs Vector containing chromosome segments
#
#  @return plot showing Battenberg average copy number profile using base plotting function
#
# create.bb.plot.average = function(BAFvals, subclones, bafsegmented, ploidy, rho, goodnessOfFit, pos_min, pos_max, segment_states_min, segment_states_tot, chr.segs) {
#   maintitle = paste("Ploidy: ", sprintf("%1.2f", ploidy), ", purity: ", sprintf("%2.0f", rho*100), "%, goodness of fit: ", sprintf("%2.1f", goodnessOfFit*100), "%", sep="")
#   nTotal = array(NA, nrow(BAFvals))
#   nMinor = array(NA, nrow(BAFvals))
#   for (i in 1:nrow(subclones)) {
#     segm_chr = subclones$chr[i] == BAFvals$Chromosome & subclones$startpos[i] < BAFvals$Position & subclones$endpos[i] >= BAFvals$Position
#     pos_min = min(which(segm_chr))
#     pos_max = max(which(segm_chr))
#     nTotal = c(nTotal, segment_states_tot[pos_min[i]:pos_max[i]])
#     nMinor = c(nMinor, segment_states_min[pos_min[i]:pos_max[i]])
#   }
#   colourTotal = "#E69F00"
#   colourMinor = "#2f4f4f"
#   base.gw.plot(bafsegmented, nTotal, nMinor, colourTotal, colourMinor, maintitle, chr.segs, y_limit, textFlag)
# }

#' @title base.gw.plot
#' @description Basis for the genome-wide plots
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param nAfullPlot Total segment copy number
#' @param nBfullPlot Segment copy number minor allele
#' @param colourTotal Colour to plot total copy number
#' @param colourMinor Colour to plot minor allele
#' @param maintitle Title comprising ploidy, rho, goodness of fit
#' @param chr.segs Vector comprising chromosome segments
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param chr.names Vector giving the names of the chromosomes as displayed on the figure
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#' @param twoColours Optional flag to specify colours, if TRUE colour is paler for CN values > y_limit
#' @keywords internal
#' @return basic plot containing chromosome positions and names, plots copy number for either ASCAT non rounded or BB average
#' @export

base.gw.plot = function(bafsegmented, nAfullPlot, nBfullPlot, colourTotal, colourMinor, maintitle, chr.segs, lrr, chr.names, y_limit, twoColours=FALSE) {
  par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  ticks=seq(0, y_limit, 1)
  plot(c(1, length(nAfullPlot)), c(0, y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
  axis(side = 2, at = ticks)
  abline(h=ticks, col="lightgrey", lty=1)

  A_rle<-rle(nAfullPlot)
  start=0
  #plot total copy number
  for (i in 1:length(A_rle$values)) {
    val<-A_rle$values[i]
    size<-A_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal, red.f=0.75, green.f=0.75, blue.f=0.75), colourTotal), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal, red.f=0.75, green.f=0.75, blue.f=0.75), colourTotal))
    start=start+size
  }

  B_rle<-rle(nBfullPlot)
  start=0
  #plot minor allele copy number
  for (i in 1:length(B_rle$values)) {
    val<-B_rle$values[i]
    size<-B_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75, green.f=0.75, blue.f=0.75), colourMinor), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75, green.f=0.75, blue.f=0.75), colourMinor))
    start=start+size
  }

  chrk_tot_len = 0
  abline(v=0, lty=1, col="lightgrey")
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]]
    chrk_hetero = intersect(names(lrr)[chrk], names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2
    text(tpos, y_limit, chr.names[i], pos = 1, cex = 2)
    abline(v=vpos, lty=1, col="lightgrey")
  }

  #   #add text to too high fragments
  #   if (textFlag) {
  #     #      rleB<-rle(nBfullPlot>y_limit)
  #     #      pos<-0
  #     #      for (i in 1:length(rleB$values)) {
  #     #        if (rleB$values[i]) {
  #     #          xpos=pos+(rleB$lengths[i]/2)
  #     #          text(xpos, y_limit+0.1, sprintf("%1.2f", nBfull[pos+1]), pos = 1, cex = 0.7)
  #     #        }
  #     #        pos=pos+rleB$lengths[i]
  #     #      }
  #
  #     rleA<-rle(nAfullPlot>y_limit)
  #     pos<-0
  #     for (i in 1:length(rleA$values)) {
  #       if (rleA$values[i]) {
  #         xpos=pos+(rleA$lengths[i]/2)
  #         text(xpos, y_limit+0.1, sprintf("%1.2f", (nAfull[pos+1]+nBfull[pos+1])), pos = 1, cex = 0.7)
  #       }
  #       pos=pos+rleA$lengths[i]
  #     }
  #   }
}

#' @title ascat.plotAscatProfile
#'
#' @description Function plotting the rounded ASCAT profiles over all chromosomes
#' @param n1all copy number major allele
#' @param n2all copy number minor allele
#' @param heteroprobes probes with heterozygous germline
#' @param ploidy ploidy of the sample
#' @param rho purity of the sample
#' @param goodnessOfFit estimated goodness of fit
#' @param nonaberrant boolean flag denoting non-aberrated samples
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#' @param ch a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22, "X"))
#'
#' @return plot showing the ASCAT profile of the sample
#'
#' @export
#'
ascat.plotAscatProfile<-function(n1all, n2all, heteroprobes, ploidy, rho, goodnessOfFit, nonaberrant, y_limit=5, ch, lrr, bafsegmented, chrs) {
  nA2 = n1all[heteroprobes]
  nB2 = n2all[heteroprobes]
  nA = ifelse(nA2>nB2, nA2, nB2)
  nB = ifelse(nA2>nB2, nB2, nA2)

  nBPlot<-ifelse(nB<=y_limit, nB+0.1, y_limit+0.1)
  nAPlot<-ifelse(nA<=y_limit, nA-0.1, y_limit+0.1)

  colourTotal="#E03546" # red
  colourMinor="#3557E0" # blue

  maintitle = paste("Ploidy: ", sprintf("%1.2f", ploidy), ", purity: ", sprintf("%2.0f", rho*100), "%, goodness of fit: ", sprintf("%2.1f", goodnessOfFit), "%", ifelse(nonaberrant, ", non-aberrant", ""), sep="")
  base.gw.plot(bafsegmented, nAPlot, nBPlot, colourTotal, colourMinor, maintitle, ch, lrr, chrs, y_limit, twoColours=TRUE)
}

#' ascat.plotGenotypes
#'
#' @param ASCATobj an ASCAT object
#' @param title main title of the plot
#' @param Tumor_BAF_noNA B-allele frequencies of the tumour sample with removed NA values
#' @param Hom Boolean vector denoting homozygous SNPs
#' @param ch_noNA vector of probes per chromosome (NA values excluded)
#'
#' @return plot showing classified BAF per sample, with unused SNPs in green, germline homozygous SNPs in blue and all others in red
#' @export
#'
ascat.plotGenotypes<-function(ASCATobj, title, Tumor_BAF_noNA, Hom, ch_noNA) {
  par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
  plot(c(1, length(Tumor_BAF_noNA)), c(0, 1), type = "n", xaxt = "n", main = title, xlab = "", ylab = "")
  points(Tumor_BAF_noNA, col=ifelse(is.na(Hom), "green", ifelse(Hom, "blue", "red")))
  points(Tumor_BAF_noNA, col=ifelse(is.na(Hom), "#00990011", ifelse(Hom, "#00005511", "#77000011")))

  abline(v=0.5, lty=1, col="lightgrey")
  chrk_tot_len = 0
  for (j in 1:length(ch_noNA)) {
    chrk = ch_noNA[[j]]
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk)
    vpos = chrk_tot_len
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2
    text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
    abline(v=vpos+0.5, lty=1, col="lightgrey")
  }
}

#' @title ascat.plotAdjustedAscatProfile
#'
#' @description Function plotting the "adjusted" (with realistic chromosome sizes) rounded/unrounded ASCAT profiles over all chromosomes.
#' @param ASCAT_output_object R object generated by the ascat.runAscat function.
#' @param REF Can be either "hg19" or "hg38" for standard human genome or a data.frame with three columns: chrom, start and end.
#' @param y_limit Optional parameter determining the size of the y axis in the profile (default=5).
#' @param plot_unrounded Optional parameter to define whether rounded (default) or unrounded profile (set to TRUE) should be plotted.
#' @param png_prefix Optional parameter to add a prefix to png name (can be also used to set a path).
#'
#' @return Plot showing the adjusted (rounded/unrounded) ASCAT profile of the sample
#'
#' @export
#'
ascat.plotAdjustedAscatProfile=function(ASCAT_output_object, REF, y_limit=5, plot_unrounded=FALSE, png_prefix="") {
  if (plot_unrounded) {
    SEGMENTS=ASCAT_output_object$segments_raw[, c(1:4, 7:8)]
    colnames(SEGMENTS)[5:6]=c("nMajor", "nMinor")
    SEGMENTS$nMajor=SEGMENTS$nMajor+SEGMENTS$nMinor
    colourA = "#943CC3" # purple
    colourB = "#60AF36" # green
  } else {
    SEGMENTS=ASCAT_output_object$segments
    SEGMENTS$nMajor=SEGMENTS$nMajor-0.1
    SEGMENTS$nMinor=SEGMENTS$nMinor+0.1
    colourA = "#E03546" # red
    colourB = "#3557E0" # blue
  }
  SEGMENTS$nMajor=ifelse(SEGMENTS$nMajor>y_limit, y_limit+0.1, SEGMENTS$nMajor)
  SEGMENTS$nMinor=ifelse(SEGMENTS$nMinor>y_limit, y_limit+0.1, SEGMENTS$nMinor)

  if (REF=="hg19") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                         135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566, 155270560))
  } else if (REF=="hg38") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                         133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                         58617616, 64444167, 46709983, 50818468, 156040895))
  } else {
    stopifnot(is.data.frame(REF))
    stopifnot(identical(colnames(REF), c("chrom", "start", "end")))
  }

  SEGMENTS$chr=gsub("^chr", "", SEGMENTS$chr)
  stopifnot(all(ASCAT_output_object$segments$chr %in% REF$chrom))
  REF$size=REF$end-REF$start+1
  REF$middle=0
  for (i in 1:nrow(REF)) {
    if (i==1) {
      REF$middle[i]=REF$size[i]/2
    } else {
      REF$middle[i]=sum(as.numeric(REF$size[1:(i-1)]))+REF$size[i]/2
    }
  }; rm(i)
  REF$cumul=cumsum(as.numeric(REF$size))
  REF$add=cumsum(as.numeric(c(0, REF$size[1:(nrow(REF)-1)])))

  SEGMENTS$startpos_adjusted=SEGMENTS$startpos
  SEGMENTS$endpos_adjusted=SEGMENTS$endpos
  for (CHR in unique(REF$chrom)) {
    INDEX=which(SEGMENTS$chr==CHR)
    if (length(INDEX)>0) {
      SEGMENTS$startpos_adjusted[INDEX]=SEGMENTS$startpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
      SEGMENTS$endpos_adjusted[INDEX]=SEGMENTS$endpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
    }
    rm(INDEX)
  }; rm(CHR)

  for (SAMPLE in sort(unique(SEGMENTS$sample))) {
    SEGS=SEGMENTS[which(SEGMENTS$sample==SAMPLE), ]
    if (nrow(SEGS)==0) warning(paste0("No segments for sample: ", SAMPLE))
    maintitle = paste("Ploidy: ", sprintf("%1.2f", ASCAT_output_object$ploidy[SAMPLE]), ", purity: ", sprintf("%2.0f", ASCAT_output_object$purity[SAMPLE]*100), "%, goodness of fit: ", sprintf("%2.1f", ASCAT_output_object$goodnessOfFit[SAMPLE]), "%", ifelse(isTRUE(ASCAT_output_object$nonaberrantarrays[SAMPLE]), ", non-aberrant", ""), sep="")
    png(filename = paste0(png_prefix, SAMPLE, ".adjusted", ifelse(plot_unrounded, "rawprofile", "ASCATprofile"), ".png"), width = 2000, height = (y_limit*100), res = 200)
    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
    ticks=seq(0, y_limit, 1)
    plot(c(1, REF$cumul[nrow(REF)]), c(0, y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
    axis(side = 2, at = ticks)
    abline(h=ticks, col="lightgrey", lty=1)
    rect(SEGS$startpos_adjusted, (SEGS$nMajor-0.07), SEGS$endpos_adjusted, (SEGS$nMajor+0.07), col=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA), border=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA))
    rect(SEGS$startpos_adjusted, (SEGS$nMinor-0.07), SEGS$endpos_adjusted, (SEGS$nMinor+0.07), col=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB), border=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB))
    abline(v=c(1, REF$cumul), lty=1, col="lightgrey")
    text(REF$middle, y_limit, REF$chrom, pos = 1, cex = 2)
    dev.off()
    rm(SEGS, ticks, maintitle)
  }; rm(SAMPLE)
}
