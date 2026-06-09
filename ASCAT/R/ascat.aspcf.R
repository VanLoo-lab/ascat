#' @title ascat.aspcf
#' @description run ASPCF segmentation
#' @details This function can be easily parallelised by controlling the selectsamples parameter\cr
#' it saves the results in LogR_PCFed[sample]_[segment].txt and BAF_PCFed[sample]_[segment].txt
#' @param ASCATobj an ASCAT object
#' @param selectsamples a vector containing the sample number(s) to PCF. Default = all
#' @param ascat.gg germline genotypes (NULL if germline data is available)
#' @param penalty penalty of introducing an additional ASPCF breakpoint (expert parameter, don't adapt unless you know what you're doing)
#' @param out.dir directory in which output files will be written. Can be set to NA to not write PCFed files.
#' @param out.prefix prefix for output file names
#' @param seed A seed to be set when subsampling SNPs for X in males (optional, default=as.integer(Sys.time())).
#'
#' @return output: ascat data structure containing:\cr
#' 1. Tumor_LogR data matrix\cr
#' 2. Tumor_BAF data matrix\cr
#' 3. Tumor_LogR_segmented: matrix of LogR segmented values\cr
#' 4. Tumor_BAF_segmented: list of BAF segmented values; each element in the list is a matrix containing the segmented values for one sample (only for probes that are not germline homozygous)\cr
#' 5. Germline_LogR data matrix\cr
#' 6. Germline_BAF data matrix\cr
#' 7. SNPpos: position of all SNPs\cr
#' 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]], ] will output the Tumor_LogR data of chromosome 13\cr
#' 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)\cr
#'
#' @export
#'
ascat.aspcf = function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL, penalty = 70, out.dir=".", out.prefix="", seed=as.integer(Sys.time())) {
  set.seed(seed)
  # first, set germline genotypes
  gg = NULL
  if (!is.null(ascat.gg)) {
    gg = ascat.gg$germlinegenotypes
  } else {
    gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 0.7
  }
  # calculate germline homozygous stretches for later resegmentation
  ghs = predictGermlineHomozygousStretches(ASCATobj$chr, gg)

  segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140))
  segmentlengths = segmentlengths[segmentlengths>=penalty]

  Tumor_LogR_segmented = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  rownames(Tumor_LogR_segmented) = rownames(ASCATobj$Tumor_LogR)
  colnames(Tumor_LogR_segmented) = colnames(ASCATobj$Tumor_LogR)
  Tumor_BAF_segmented = list()
  for (sample in selectsamples) {
    print.noquote(paste("Sample ", ASCATobj$samples[sample], " (", sample, "/", length(ASCATobj$samples), ")", sep=""))
    logrfilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], ".LogR.PCFed.txt", sep=""))
    baffilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], ".BAF.PCFed.txt", sep=""))
    logRPCFed = numeric(0)
    bafPCFed = numeric(0)
    # specific process for nonPAR in males
    if (!is.null(ASCATobj$X_nonPAR) && ASCATobj$gender[sample]=="XY") {
      # select SNPs with non-NA BAF values in nonPAR region
      nonPAR_index=which(ASCATobj$SNPpos[, 1] %in% c("chrX", "X") & ASCATobj$SNPpos[, 2]>=ASCATobj$X_nonPAR[1] & ASCATobj$SNPpos[, 2]<=ASCATobj$X_nonPAR[2] & !is.na(gg[, sample]))
      # store hmz/htz information for autosomes
      autosomes_info=table(gg[which(ASCATobj$SNPpos[, 1] %in% setdiff(ASCATobj$chrs, ASCATobj$sexchromosomes)), sample])
      if (length(nonPAR_index)>5) {
        # set all to hmz
        gg[nonPAR_index, sample]=TRUE
        if (!is.null(ASCATobj$Germline_BAF)) {
          # compute distance to BAF=0/1
          DIST=1-sapply(ASCATobj$Germline_BAF[nonPAR_index, sample], function(x) {if (x>0.5) return(x) else return(1-x)})
          # select X% (derived from autosomes) of SNPs based on closest distance to BAF=0/1 and force those to be considered for ASPCF
          gg[nonPAR_index[which(rank(DIST, ties.method="random")<=round(length(DIST) * (autosomes_info["FALSE"]/sum(autosomes_info))))], sample]=FALSE
          rm(DIST)
        } else {
          gg[sample(nonPAR_index, round(length(nonPAR_index) * (autosomes_info["FALSE"]/sum(autosomes_info)))), sample]=FALSE
        }
      }
      rm(nonPAR_index, autosomes_info)
    }
    for (segmentlength in segmentlengths) {
      logRPCFed = numeric(0)
      bafPCFed = numeric(0)
      tbsam = ASCATobj$Tumor_BAF[, sample]
      names(tbsam) = rownames(ASCATobj$Tumor_BAF)
      homosam = gg[, sample]
      for (chrke in 1:length(ASCATobj$chr)) {
        lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]], sample]
        #winsorize to remove outliers
        #this has a problem with NAs
        lrwins = vector(mode="numeric", length=length(lr))
        lrwins[is.na(lr)] = NA
        lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)], 2.5, 25)$ywin
        baf = tbsam[ASCATobj$chr[[chrke]]]
        homo = homosam[ASCATobj$chr[[chrke]]]
        Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
        bafsel = baf[Select_het]
        # winsorize BAF as well (as a safeguard)
        bafselwinsmirrored = madWins(ifelse(bafsel>0.5, bafsel, 1-bafsel), 2.5, 25)$ywin
        bafselwins = ifelse(bafsel>0.5, bafselwinsmirrored, 1-bafselwinsmirrored)
        indices = which(Select_het)
        logRaveraged = NULL
        if (length(indices)!=0) {
          averageIndices = c(1, (indices[1:(length(indices)-1)]+indices[2:length(indices)])/2, length(lr)+0.01)
          startindices = ceiling(averageIndices[1:(length(averageIndices)-1)])
          endindices = floor(averageIndices[2:length(averageIndices)]-0.01)
          if (length(indices)==1) {
            startindices = 1
            endindices = length(lr)
          }
          #nrIndices = endindices - startindices + 1
          logRaveraged = vector(mode="numeric", length=length(indices))
          for (i in 1:length(indices)) {
            if (is.na(endindices[i])) {
              endindices[i]=startindices[i]
            }
            logRaveraged[i]=mean(lrwins[startindices[i]:endindices[i]], na.rm=TRUE)
          }
        }
        # if there are no probes in the segment (after germline homozygous removal), don't do anything, except add a LogR segment
        if (length(logRaveraged)>0) {
          logRASPCF = NULL
          bafASPCF = NULL
          if (length(logRaveraged)<6) {
            logRASPCF = rep(mean(logRaveraged), length(logRaveraged))
            bafASPCF = rep(mean(bafselwinsmirrored), length(logRaveraged))
            if ("isTargetedSeq" %in% names(ASCATobj) && ASCATobj$isTargetedSeq && bafASPCF[1]<=0.55) bafASPCF=rep(0.5, length(logRaveraged))
          } else {
            PCFed = fastAspcf(logRaveraged, bafselwins, 6, segmentlength, ASCATobj$isTargetedSeq)
            logRASPCF = PCFed$yhat1
            bafASPCF = PCFed$yhat2
          }
          names(bafASPCF)=names(indices)
          logRc = numeric(0)
          for (probe in 1:length(logRASPCF)) {
            if (probe == 1) {
              logRc = rep(logRASPCF[probe], indices[probe])
            }
            # if probe is 1, set the beginning, and let the loop go:
            if (probe == length(logRASPCF)) {
              logRc = c(logRc, rep(logRASPCF[probe], length(lr)-indices[probe]))
            } else if (logRASPCF[probe]==logRASPCF[probe+1]) {
              logRc = c(logRc, rep(logRASPCF[probe], indices[probe+1]-indices[probe]))
            } else {
              #find best breakpoint
              d = numeric(0)
              totall = indices[probe+1]-indices[probe]
              for (bp in 0:(totall-1)) {
                dis = sum(abs(lr[(1:bp)+indices[probe]]-logRASPCF[probe]), na.rm=TRUE)
                if (bp!=totall) {
                  dis = sum(dis, sum(abs(lr[((bp+1):totall)+indices[probe]]-logRASPCF[probe+1]), na.rm=TRUE), na.rm=TRUE)
                }
                d = c(d, dis)
              }
              breakpoint = which.min(d)-1
              logRc = c(logRc, rep(logRASPCF[probe], breakpoint), rep(logRASPCF[probe+1], totall-breakpoint))
            }
          }
          #2nd step: adapt levels!
          logRd = numeric(0)
          seg = rle(logRc)$lengths
          startprobe = 1
          endprobe = 0
          for (i in 1:length(seg)) {
            endprobe = endprobe+seg[i]
            level = mean(lr[startprobe:endprobe], na.rm=TRUE)
            logRd = c(logRd, rep(level, seg[i]))
            startprobe = startprobe + seg[i]
          }
          logRPCFed = c(logRPCFed, logRd)
          bafPCFed = c(bafPCFed, bafASPCF)
        } else { # add a LogR segment
          level = mean(lr, na.rm=TRUE)
          reps = length(lr)
          logRPCFed = c(logRPCFed, rep(level, reps))
        }
        # correct wrong segments in germline homozygous stretches:
        homsegs = ghs[[sample]][ghs[[sample]][, 1]==chrke, ]
        startchr = min(ASCATobj$chr[[chrke]])
        endchr = max(ASCATobj$chr[[chrke]])
        # to solve an annoying error when homsegs has length 1:
        if (length(homsegs)==3) {
          homsegs=t(as.matrix(homsegs))
        }
        if (!is.null(homsegs)&&dim(homsegs)[1]!=0) {
          for (i in 1:dim(homsegs)[1]) {
            if (anyNA(homsegs[i, ])) next
            # note that only the germline homozygous segment is resegmented, plus a bit extra (but this is NOT replaced)
            #startpos = max(homsegs[i, 2], startchr)
            #endpos = min(homsegs[i, 3], endchr)
            # PCF over a larger fragment
            startpos2 = max(homsegs[i, 2]-100, startchr)
            endpos2 = min(homsegs[i, 3]+100, endchr)
            # take into account a little extra (difference between startpos2 and startpos3 is not changed)
            startpos3 = max(homsegs[i, 2]-5, startchr)
            endpos3 = min(homsegs[i, 3]+5, endchr)
            # note that the parameters are arbitrary, but <100 seems to work on the ERBB2 example!
            # segmentlength is lower here, as in the full data, noise on LogR is higher!
            # do this on Winsorized data too!
            towins = ASCATobj$Tumor_LogR[startpos2:endpos2, sample]
            winsed = madWins(towins[!is.na(towins)], 2.5, 25)$ywin
            pcfed = vector(mode="numeric", length=length(towins))
            pcfed[!is.na(towins)] = exactPcf(winsed, 6, floor(segmentlength/4))
            pcfed2 = pcfed[(startpos3-startpos2+1):(endpos3-startpos2+1)]
            dif = abs(pcfed2-logRPCFed[startpos3:endpos3])
            #0.3 is hardcoded here, in order not to have too many segments!
            #only replace if enough probes differ (in order not to get singular probes with different breakpoints)
            if (!anyNA(dif)&&sum(dif>0.3)>5) {
              #replace a bit more to make sure no 'lone' probes are left (startpos3 instead of startpos)
              logRPCFed[startpos3:endpos3]=ifelse(dif>0.3, pcfed2, logRPCFed[startpos3:endpos3])
            }
          }
        }
      }
      #fill in NAs (otherwise they cause problems):
      #some NA probes are filled in with zero, replace those too:
      logRPCFed = fillNA(logRPCFed, zeroIsNA=TRUE)

      #adapt levels again
      seg = rle(logRPCFed)$lengths
      logRPCFed = numeric(0)
      startprobe = 1
      endprobe = 0
      prevlevel = 0
      for (i in 1:length(seg)) {
        endprobe = endprobe+seg[i]
        level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe, sample], na.rm=TRUE)
        #making sure no NA's get filled in...
        if (is.nan(level)) {
          level=prevlevel
        } else {
          prevlevel=level
        }
        logRPCFed = c(logRPCFed, rep(level, seg[i]))
        startprobe = startprobe + seg[i]
      }
      #put in names and write results to files
      names(logRPCFed) = rownames(ASCATobj$Tumor_LogR)

      # if less than 800 segments: this segmentlength is ok, otherwise, rerun with higher segmentlength
      if (length(unique(logRPCFed))<800) {
        break
      }
    }

    if (!is.na(out.dir)) write.table(logRPCFed, logrfilename, sep="\t", col.names=FALSE)
    if (!is.na(out.dir)) write.table(bafPCFed, baffilename, sep="\t", col.names=FALSE)
    bafPCFed = as.matrix(bafPCFed)
    Tumor_LogR_segmented[, sample] = logRPCFed
    Tumor_BAF_segmented[[sample]] = 1-bafPCFed
  }

  ASCATobj$Tumor_LogR_segmented=Tumor_LogR_segmented
  ASCATobj$Tumor_BAF_segmented=Tumor_BAF_segmented
  ASCATobj$failedarrays=ascat.gg$failedarrays
  return(ASCATobj)
}

#' @title predictGermlineHomozygousStretches
#' @description helper function to predict germline homozyguous segments for later re-segmentation
#' @param chr a list containing vectors with the indices for each distinct part that can be segmented separately
#' @param hom germline genotypes
#'
#' @keywords internal
#'
#' @return germline homozyguous segments
#'
#'
predictGermlineHomozygousStretches = function(chr, hom) {

  # contains the result: a list of vectors of probe numbers in homozygous stretches for each sample
  HomoStretches = list()

  for (sam in 1:dim(hom)[2]) {
    homsam = hom[, sam]

    perchom = sum(homsam, na.rm=TRUE)/sum(!is.na(homsam))

    # NOTE THAT A P-VALUE THRESHOLD OF 0.001 IS HARDCODED HERE
    homthres = ceiling(log(0.001, perchom))

    allhprobes = NULL
    for (chrke in 1:length(chr)) {
      hschr = homsam[chr[[chrke]]]

      hprobes = vector(length=0)
      for (probe in 1:length(hschr)) {
        if (!is.na(hschr[probe])) {
          if (hschr[probe]) {
            hprobes = c(hprobes, probe)
          } else {
            if (length(hprobes)>=homthres) {
              allhprobes = rbind(allhprobes, c(chrke, chr[[chrke]][min(hprobes)], chr[[chrke]][max(hprobes)]))
            }
            hprobes = vector(length=0)
          }
        }
      }
      # if the last probe is homozygous, this is not yet accounted for
      if (!is.na(hschr[probe]) && hschr[probe]) {
        if (length(hprobes)>=homthres) {
          allhprobes = rbind(allhprobes, c(chrke, chr[[chrke]][min(hprobes)], chr[[chrke]][max(hprobes)]))
        }
      }

    }

    if (is.null(allhprobes)) allhprobes=rbind(NULL, c(0, 0, 0))
    HomoStretches[[sam]]=allhprobes

  }

  return(HomoStretches)
}

#
# Enhanced bivariate PCF filter for aCGH data (v. 08.02.2010)
# Whole chromosomes/chromosome arms wrapper function
#

fastAspcf <- function(logR, allB, kmin, gamma, isTargetedSeq) {
  if (is.null(isTargetedSeq)) isTargetedSeq=FALSE

  N <- length(logR)
  w <- 1000 #w: windowsize
  d <- 100

  startw = -d
  stopw = w-d

  nseg = 0
  var2 = 0
  var3 = 0
  breakpts = 0
  larger = TRUE
  repeat {
    from <- max(c(1, startw))
    to  <-  min(c(stopw, N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]

    sd1 <- getMad(logRpart)
    sd2 <- getMad(allBflip)
    sd3 <- getMad(allBpart)

    #Must check that sd1 and sd2 are defined and != 0:
    sd.valid <- c(!is.na(sd1), !is.na(sd2), sd1!=0, sd2!=0)
    if (all(sd.valid)) {
      #run aspcfpart:
      #part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      breakptspart <- part.res$breakpts
      # the 'larger' is (occasionally) necessary in the last window of the segmentation!
      larger = breakptspart>breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      var3 <- var3 + sd3^2
      nseg = nseg+1
    }

    if (stopw < N+d) {
      startw <- min(stopw-2*d + 1, N-2*d)
      stopw <- startw + w
    } else {
      break
    }

  }#end repeat
  breakpts <- unique(c(breakpts, N))
  if (nseg==0) { nseg=1 }  #just in case the sd-test never passes.
  sd2 <- sqrt(var2/nseg)
  sd3 <- sqrt(var3/nseg)

  # On each segment calculate mean of unflipped B allele data
  frst <- breakpts[1:length(breakpts)-1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)

  yhat1 <- rep(NA, N)
  yhat2 <- rep(NA, N)

  for (i in 1:nseg) {
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i]-frst[i]+1)
    yi2 <- allB[frst[i]:last[i]]
    # Center data around zero (by subtracting 0.5) and estimate mean
    if (length(yi2)== 0) {
      mu <- 0
    } else {
      mu <- mean(abs(yi2-0.5))
    }

    # Make a (slightly arbitrary) decision concerning branches
    # This may be improved by a test of equal variances
    if (sqrt(sd2^2+mu^2) < 2*sd2) {
      # if (sd3 < 1.8*sd2) {
      mu <- 0
    }
    if (isTargetedSeq && mu<=0.05) mu=0
    yhat2[frst[i]:last[i]] <- rep(mu+0.5, last[i]-frst[i]+1)
  }
  return(list(yhat1=yhat1, yhat2=yhat2))

}#end fastAspcf



aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma) {

  from <- max(c(1, a))
  usefrom <- max(c(1, a+d))
  useto <- min(c(N, b-d))

  N <- length(logRpart)
  y1 <- logRpart
  y2 <- allBflip

  #Check that vectors are long enough to run algorithm:
  if (N < 2*kmin) {
    breakpts <- 0
    return(list(breakpts=breakpts))
  }

  # Find initSum, initKvad, initAve for segment y[1..kmin]
  initSum1 <- sum(y1[1:kmin])
  initKvad1 <- sum(y1[1:kmin]^2)
  initAve1 <- initSum1/kmin
  initSum2 <- sum(y2[1:kmin])
  initKvad2 <- sum(y2[1:kmin]^2)
  initAve2 <- initSum2/kmin

  # Define vector of best costs
  bestCost <- rep(0, N)
  cost1 <- (initKvad1 - initSum1*initAve1)/sd1^2
  cost2 <- (initKvad2 - initSum2*initAve2)/sd2^2
  bestCost[kmin] <- cost1 + cost2

  # Define vector of best splits
  bestSplit <- rep(0, N)

  # Define vector of best averages
  bestAver1 <- rep(0, N)
  bestAver2 <- rep(0, N)
  bestAver1[kmin] <- initAve1
  bestAver2[kmin] <- initAve2


  #Initialize
  Sum1 <- rep(0, N)
  Sum2 <- rep(0, N)
  Kvad1 <- rep(0, N)
  Kvad2 <- rep(0, N)
  Aver1 <- rep(0, N)
  Aver2 <- rep(0, N)
  Cost <- rep(0, N)

  # We have to treat the region y(1..2*kmin-1) separately, as it
  # cannot be split into two full segments
  kminP1 <- kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum1[kminP1:k] <- Sum1[kminP1:k]+y1[k]
    Aver1[kminP1:k] <- Sum1[kminP1:k] / ((k-kmin):1)
    Kvad1[kminP1:k] <- Kvad1[kminP1:k]+y1[k]^2
    Sum2[kminP1:k] <- Sum2[kminP1:k]+y2[k]
    Aver2[kminP1:k] <- Sum2[kminP1:k] / ((k-kmin):1)
    Kvad2[kminP1:k] <- Kvad2[kminP1:k]+y2[k]^2


    bestAver1[k] <- (initSum1+Sum1[kminP1])/k
    bestAver2[k] <- (initSum2+Sum2[kminP1])/k
    cost1 <- ((initKvad1+Kvad1[kminP1])-k*bestAver1[k]^2)/sd1^2
    cost2 <- ((initKvad2+Kvad2[kminP1])-k*bestAver2[k]^2)/sd2^2

    bestCost[k] <- cost1 + cost2

  }


  for (n in (2*kmin):N) {

    nMkminP1=n-kmin+1

    Sum1[kminP1:n] <- Sum1[kminP1:n]+ y1[n]
    Aver1[kminP1:n] <- Sum1[kminP1:n] / ((n-kmin):1)
    Kvad1[kminP1:n] <- Kvad1[kminP1:n]+ (y1[n])^2

    cost1 <- (Kvad1[kminP1:nMkminP1]-Sum1[kminP1:nMkminP1]*Aver1[kminP1:nMkminP1])/sd1^2

    Sum2[kminP1:n] <- Sum2[kminP1:n]+ y2[n]
    Aver2[kminP1:n] <- Sum2[kminP1:n] / ((n-kmin):1)
    Kvad2[kminP1:n] <- Kvad2[kminP1:n]+ (y2[n])^2
    cost2 <- (Kvad2[kminP1:nMkminP1]-Sum2[kminP1:nMkminP1]*Aver2[kminP1:nMkminP1])/sd2^2

    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2

    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos] + gamma

    aver1 <- Aver1[Pos]
    aver2 <- Aver2[Pos]
    totAver1 <- (Sum1[kminP1]+initSum1)/n
    totCost1 <- ((Kvad1[kminP1]+initKvad1) - n*totAver1*totAver1)/sd1^2
    totAver2 <- (Sum2[kminP1]+initSum2)/n
    totCost2 <- ((Kvad2[kminP1]+initKvad2) - n*totAver2*totAver2)/sd2^2
    totCost <- totCost1 + totCost2

    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver1 <- totAver1
      aver2 <- totAver2
    }
    bestCost[n] <- cost
    bestAver1[n] <- aver1
    bestAver2[n] <- aver2
    bestSplit[n] <- Pos-1


  }#endfor


  # Trace back
  n <- N
  breakpts <- n
  while (n > 0) {
    breakpts <- c(bestSplit[n], breakpts)
    n <- bestSplit[n]
  }#endwhile

  breakpts <- breakpts + from -1
  breakpts <- breakpts[breakpts>=usefrom & breakpts<=useto]

  return(list(breakpts=breakpts))

}#end aspcfpart

#Get mad SD (based on KL code)
getMad <- function(x, k=25) {

  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]

  #Calculate runMedian
  runMedian <- medianFilter(x, k)

  dif <- x-runMedian
  SD <- mad(dif)

  return(SD)
}

exactPcf <- function(y, kmin, gamma) {
  ## Implementaion of exact PCF by Potts-filtering
  ## x: input array of (log2) copy numbers
  ## kmin: Mininal length of plateaus
  ## gamma: penalty for each discontinuity
  N <- length(y)
  yhat <- rep(0, N)
  if (N < 2*kmin) {
    yhat <- rep(mean(y), N)
    return(yhat)
  }
  initSum <- sum(y[1:kmin])
  initKvad <- sum(y[1:kmin]^2)
  initAve <- initSum/kmin
  bestCost <- rep(0, N)
  bestCost[kmin] <- initKvad - initSum*initAve
  bestSplit <- rep(0, N)
  bestAver <- rep(0, N)
  bestAver[kmin] <- initAve
  Sum <- rep(0, N)
  Kvad <- rep(0, N)
  Aver <- rep(0, N)
  Cost <- rep(0, N)
  kminP1=kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    Aver[kminP1:k] <- Sum[kminP1:k] / ((k-kmin):1)
    Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
    bestAver[k] <- (initSum+Sum[kminP1])/k
    bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
  }
  for (n in (2*kmin):N) {
    yn <- y[n]
    yn2 <- yn^2
    Sum[kminP1:n] <- Sum[kminP1:n]+yn
    Aver[kminP1:n] <- Sum[kminP1:n] / ((n-kmin):1)
    Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
    nMkminP1=n-kmin+1
    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos]
    aver <- Aver[Pos]
    totAver <- (Sum[kminP1]+initSum)/n
    totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver <- totAver
    }
    bestCost[n] <- cost
    bestAver[n] <- aver
    bestSplit[n] <- Pos-1
  }
  n <- N
  while (n > 0) {
    yhat[(bestSplit[n]+1):n] <- bestAver[n]
    n <- bestSplit[n]
  }
  return(yhat)
}


#Perform MAD winsorization:
madWins <- function(x, tau, k) {
  xhat <- medianFilter(x, k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin, sdev=SD, outliers=outliers))
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x, k) {
  n <- length(x)
  filtWidth <- 2*k + 1

  #Make sure filtWidth does not exceed n
  if (filtWidth > n) {
    if (n==0) {
      filtWidth <- 1
    }else if (n %% 2 == 0) {
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    } else {
      filtWidth <- n
    }
  }

  runMedian <- runmed(x, k=filtWidth, endrule="median")

  return(runMedian)

}

fillNA = function(vec, zeroIsNA=TRUE) {
  if (zeroIsNA) {vec[vec==0] <- NA}
  nas = which(is.na(vec))

  if (length(nas) == 0) {
    return(vec)
  }

  # Find stretches of contiguous NAs
  starts = c(1, which(diff(nas)>1)+1)
  ends = c(starts[-1] - 1, length(nas))

  starts = nas[starts]
  ends = nas[ends]

  # Special-case: vec[1] is NA
  startAt = 1
  if (starts[1]==1) {
    vec[1:ends[1]] = vec[ends[1]+1]
    startAt = 2
  }

  if (startAt > length(starts)) {
    return(vec)
  }

  # Special-case: last element in vec is NA
  endAt = length(starts)
  if (is.na(vec[length(vec)])) {
    vec[starts[endAt]:ends[endAt]] = vec[starts[endAt]-1]
    endAt = endAt-1
  }

  if (endAt < startAt) {
    return(vec)
  }

  # For each stretch of NAs, set start:midpoint to the value before,
  # and midpoint+1:end to the value after.
  for (i in startAt:endAt) {
    start = starts[i]
    end = ends[i]
    N = 1 + end-start
    if (N==1) {
      vec[start] = vec[start-1]
    } else {
      midpoint = start+ceiling(N/2)
      vec[start:midpoint] = vec[start-1]
      vec[(midpoint+1):end] = vec[end+1]
    }
  }

  return(vec)
}

psi <- function(x, z) {
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}
