#' Allele-specific segmentation of multiple samples
#'
#' @description This segmentation function should only be used if part of the breakpoints are expected to be shared
#' between samples, e.g. due to a common ancestry.
#'
#' @param ASCATobj an ASCAT object
#' @param ascat.gg germline genotypes (NULL if germline data is available)
#' @param penalty penalty of introducing an additional ASPCF breakpoint (expert parameter, don't adapt unless you know what you are doing)
#' @param out.dir directory in which output files will be written. Can be set to NA to not write PCFed files.
#' @param wsample Vector of length length(ASCATobj$samples). Can be used to assign different weights to samples, for example to account for differences in sequencing quality. (Default = NULL)
#' @param selectAlg Set to "exact" to run the exact algorithm, or "fast" to run the heuristic algorithm. (Default = "exact")
#' @param refine Logical. Should breakpoints be refined on a per sample base? Otherwise each breakpoint is assumed to be present in each sample. (Default = TRUE)
#' @param seed A seed to be set when subsampling SNPs for X in males (optional, default=as.integer(Sys.time())).
#'
#' @details This function saves the results in in [sample].LogR.PCFed.txt and [sample].BAF.PCFed.txt
#'
#' @return output: ascat data structure containing:\cr
#' 1. Tumor_LogR data matrix\cr
#' 2. Tumor_BAF data matrix\cr
#' 3. Tumor_LogR_segmented: matrix of LogR segmented values\cr
#' 4. Tumor_BAF_segmented: list of BAF segmented values; each element in the list is a matrix containing the segmented values for one sample (only for probes that are germline homozygous)\cr
#' 5. Germline_LogR data matrix\cr
#' 6. Germline_BAF data matrix\cr
#' 7. SNPpos: position of all SNPs\cr
#' 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]], ] will output the Tumor_LogR data of chromosome 13\cr
#' 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)\cr
#'
#' @export
#'
ascat.asmultipcf <- function(ASCATobj, ascat.gg = NULL, penalty = 70, out.dir = ".", wsample=NULL, selectAlg="exact", refine=TRUE, seed=as.integer(Sys.time())) {
  if (is.null(ASCATobj$isTargetedSeq)) ASCATobj$isTargetedSeq=FALSE
  set.seed(seed)
  useLogRonlySites=TRUE
  #first, set germline genotypes
  gg = NULL
  if (!is.null(ascat.gg)) {
    gg = ascat.gg$germlinegenotypes
  } else {
    gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 0.7
  }

  # in asmultipcf, we're expecting one germline (could be multiple germlines but we're only using the first one)
  if (ncol(gg)>1) gg=gg[, 1]
  # specific process for nonPAR in males
  if (!is.null(ASCATobj$X_nonPAR) && ASCATobj$gender[1]=="XY") {
    # select SNPs with non-NA BAF values in nonPAR region
    nonPAR_index=which(ASCATobj$SNPpos[, 1] %in% c("chrX", "X") & ASCATobj$SNPpos[, 2]>=ASCATobj$X_nonPAR[1] & ASCATobj$SNPpos[, 2]<=ASCATobj$X_nonPAR[2] & !is.na(gg))
    # store hmz/htz information for autosomes
    autosomes_info=table(gg[which(ASCATobj$SNPpos[, 1] %in% setdiff(ASCATobj$chrs, ASCATobj$sexchromosomes))])
    if (length(nonPAR_index)>5) {
      # set all to hmz
      gg[nonPAR_index]=TRUE
      if (!is.null(ASCATobj$Germline_BAF)) {
        # compute distance to BAF=0/1
        DIST=1-sapply(ASCATobj$Germline_BAF[nonPAR_index, 1], function(x) {if (x>0.5) return(x) else return(1-x)})
        # select X% (derived from autosomes) of SNPs based on closest distance to BAF=0/1 and force those to be considered for ASPCF
        gg[nonPAR_index[which(rank(DIST, ties.method="random")<=round(length(DIST) * (autosomes_info["FALSE"]/sum(autosomes_info))))]]=FALSE
        rm(DIST)
      } else {
        gg[sample(nonPAR_index, round(length(nonPAR_index) * (autosomes_info["FALSE"]/sum(autosomes_info))))]=FALSE
      }
    }
    rm(nonPAR_index, autosomes_info)
  }

  segmentlengths = unique(c(penalty, 25, 50, 100, 200, 400, 800))
  segmentlengths = segmentlengths[segmentlengths>=penalty]
  for (segmentlength in segmentlengths) {
    print(paste("Segmentlength", segmentlength)) ## TODO
    bafnames = character()
    logRPCFed = matrix(nrow=0, ncol=length(ASCATobj$samples))
    bafPCFed = matrix(nrow=0, ncol=length(ASCATobj$samples))
    for (chrke in 1:length(ASCATobj$chr)) {
      ## extract logR as well as tumor and normal baf of a segmentable chromosome arm (chrke)
      lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]], ]
      baf = ASCATobj$Tumor_BAF[ASCATobj$chr[[chrke]], ]
      homo = gg[ASCATobj$chr[[chrke]]]

      ## set baf of homozygous sites to missing
      bafna <- baf
      bafna[homo|is.na(homo), ] <- NA
      ## select sites used for segmentation
      if (useLogRonlySites) {
        Select_sites <- !(apply(bafna, 1, function(x) all(is.na(x))) & apply(lr, 1, function(x) all(is.na(x))))
        Select_sites2 <- apply(bafna, 1, function(x) !all(is.na(x))) & apply(lr, 1, function(x) !all(is.na(x)))
        Select_sites2 <- Select_sites2[which(Select_sites)]
      } else {
        Select_sites <- apply(bafna, 1, function(x) !all(is.na(x))) & apply(lr, 1, function(x) !all(is.na(x)))
      }

      bafsel = bafna[Select_sites, , drop=FALSE]
      indices = which(Select_sites)
      if (useLogRonlySites) {
        indices_bafOutput = which(Select_sites)
        indices_bafOutput <- indices_bafOutput[Select_sites2]
      } else {
        indices_bafOutput = which(Select_sites)
      }
      bafnames= c(bafnames, names(indices_bafOutput))
      if (length(indices)!=0) {
        ## mirror baf so that areas with allelic imbalance have only one band
        bafselmirrored <- mirrorBafMatrix(bafsel)

        ## winsorize lr and baf to remove outliers on a per sample base
        lrwins <- madWinsMatrixWithNA(lr, tau=2.5, k=25)
        bafwins <- mirrorBafMatrix(madWinsMatrixWithNA(bafselmirrored, tau=2.5, k=25))

        ## average logR values around each selected het site
        logRaveraged = NULL

        averageIndices = c(1, (indices[1:(length(indices)-1)]+indices[2:length(indices)])/2, nrow(lr)+0.01)
        startindices = ceiling(averageIndices[1:(length(averageIndices)-1)])
        endindices = floor(averageIndices[2:length(averageIndices)]-0.01)
        if (length(indices)==1) {
          startindices = 1
          endindices = nrow(lr)
        }
        #nrIndices = endindices - startindices + 1
        logRaveraged = matrix(nrow=nrow(bafwins), ncol=ncol(bafwins))

        for (i in 1:length(indices)) {
          if (is.na(endindices[i])) {
            endindices[i]=startindices[i]
          }
          logRaveraged[i, ]=apply(lrwins, 2, function(x) {
            mean(x[startindices[i]:endindices[i]], na.rm=TRUE)
          })
        }
      }

      ## SEGMENTATION
      if (length(indices)==0) {
        # if there are no probes in the segment, don't do anything, except add a LogR segment
        level = apply(lr, 2, function(x) mean(x, na.rm=TRUE))
        reps = nrow(lr)
        logRPCFed = rbind(logRPCFed, t(as.matrix(level))[rep(1, reps), ])
      } else {
        if (nrow(logRaveraged)<6) {
          if (nrow(logRaveraged)==1) {
            logRASPCF = logRaveraged
            bafASPCF = bafwins
          } else {
            logRASPCF = apply(logRaveraged, 2, function(x) rep(mean(x, na.rm=TRUE), length(x)))
            bafASPCF = apply(bafwins, 2, function(x) rep(ifelse(mean(x, na.rm=TRUE)>=0.5, mean(x, na.rm=TRUE), 1-mean(x, na.rm=TRUE)), length(x)))
            if (ASCATobj$isTargetedSeq) apply(bafASPCF, 2, function(x) ifelse(x<=0.55, 0.5, x))
          }
        } else {
          # combine logR and BAF data into one matrix for joint segmentation
          lrANDbaf <- cbind(logRaveraged, bafwins)

          # create weight matrix with same dimensions as t(lrANDbaf) where NA is 0 and everything else is 1
          w <- matrix(1, ncol=nrow(lrANDbaf), nrow=ncol(lrANDbaf))
          w[is.na(t(lrANDbaf))] = 0
          # if samples weights provided change overall weights accordingly
          if (!is.null(wsample)) {
            if (length(wsample)==ncol(lrANDbaf)) {
              w <- w*wsample
            } else if (length(wsample)==ncol(lrANDbaf)/2) {
              wsample <- c(wsample, wsample)
              w <- w*wsample
            } else {
              stop(paste("Length of parameter wsample is", length(wsample), "but should be", ncol(lrANDbaf)/2, "or", ncol(lrANDbaf)))
            }
          }
          # replace NAs in lrANDbaf with arbitrary value, e.g. 0 (so that weight * value = 0 for NA entries)
          lrANDbaf[is.na(lrANDbaf)] = 0

          # run segmentation
          if (selectAlg=="exact") {
            MULTIPCF <- ASmultiPCFcompact(nr=w, wSum=t(lrANDbaf)*w, gamma=segmentlength, wsample = wsample, nProbesPerSeg=rep(1, ncol(w))) ## takes samples in rows; wSum has to be the weighted sum!
            MULTIPCF$yhat <- expandMulti(nProbes = nrow(lrANDbaf),
                                         nSamples = ncol(lrANDbaf),
                                         lengthInt =  MULTIPCF$Lengde,
                                         mean = MULTIPCF$mean)
            names(MULTIPCF) <- c("length", "start0", "mean", "nIntervals", "yhat")
          } else if (selectAlg=="fast") {
            MULTIPCF <- runFastASMultiPCF(x=lrANDbaf, w=w, gamma=segmentlength, yest=TRUE, wsample = wsample, nProbesPerSeg=rep(1, ncol(w)))
          } else {
            stop(paste0("Unknown segmentation algorithm '", selectAlg, "'. Should be 'exact' or 'fast'."))
          }
          logRASPCF <- t(MULTIPCF$yhat)[, 1:ncol(logRaveraged)]
          bafASPCF <- t(MULTIPCF$yhat)[, (1:ncol(bafwins))+ncol(bafwins)]

          ## per-sample segmentation based on previously found multi-sample breakpoint optimization
          if (refine) {
            if (length(MULTIPCF$start0)>1) { ## segments can only be merged if there are more than 1
              for (sample in 1:length(ASCATobj$samples)) {
                ## compress LogR and BAF data for each sample based on global break points
                sampleidx <- c(sample, sample+length(ASCATobj$samples))
                compX <- compactASMulti(y = t(lrANDbaf)[sampleidx, ],
                                        mark = 1:nrow(lrANDbaf) %in% (c(MULTIPCF$start0[-1]-1, nrow(lrANDbaf))),
                                        w = w[sampleidx, ])
                ## run exact algorithm on compressed data (penalty parameter is scaled by number of samples)
                compPotts <- ASmultiPCFcompact(nr=compX$Weight, wSum=compX$wSum, gamma = segmentlength/length(ASCATobj$samples), wsample = wsample[sampleidx], nProbesPerSeg=compX$nProbesPerSeg)
                ## if segments were merged, update means
                if (compPotts$nIntervals<MULTIPCF$nIntervals) {
                  potts <- expandMulti(nProbes = nrow(lrANDbaf), nSamples = 2, lengthInt =  compPotts$Lengde, mean=compPotts$mean)
                  logRASPCF[, sample] <- potts[1, ]
                  bafASPCF[, sample] <- potts[2, ]
                }
              }
            }
          }
        }

        ## only select BAF sites whose germline genotype is heterozygous and not missing for output + not all BAF values missing across samples
        if (useLogRonlySites) {
          bafASPCF <- bafASPCF[Select_sites2, , drop=FALSE]
        }

        if (nrow(bafASPCF)>0) {
          ## correct baf estimate of balanced streches which is skewed due to mirroring (part of fastAspcf in single sample algorithm)
          bafPCFed_s = matrix(nrow=nrow(bafASPCF), ncol=ncol(bafASPCF))
          for (sample in 1:ncol(bafASPCF)) {

            ## calculate location of breakpoints
            breakpts <- c(0, cumsum(rle(bafASPCF[, sample])$lengths))

            # On each segment calculate mean of unflipped B allele data
            frst <- breakpts[1:(length(breakpts)-1)] + 1
            nseg = length(frst)
            last <- breakpts[2:length(breakpts)]

            ## get median absolute deviation from running median of baf
            sd <- getMadwithNA(bafwins[!(homo|is.na(homo))[indices], sample])

            yhat <- rep(NA, nrow(bafASPCF))

            for (i in 1:nseg) {
              yi <- bafASPCF[frst[i]:last[i], sample]

              if (!all(is.na(yi))) {

                # Center data around zero (by subtracting 0.5) and estimate mean
                if (length(yi)== 0) {
                  mu <- 0
                } else {
                  mu <- mean(abs(yi-0.5), na.rm = TRUE)
                }

                # Make a (slightly arbitrary) decision concerning branches
                # This may be improved by a test of equal variances
                if (sqrt(sd^2+mu^2) < 2*sd) {
                  mu <- 0
                }
                if (ASCATobj$isTargetedSeq && mu<=0.05) mu=0
                yhat[frst[i]:last[i]] <- rep(mu+0.5, last[i]-frst[i]+1)
              }
            }
            bafPCFed_s[, sample] = yhat
          } ## end baf correction
          bafPCFed = rbind(bafPCFed, bafPCFed_s)
        }

        ## expand logR data to full size to conform with standard package output
        ## start:
        probe=1
        logRC=logRASPCF[rep(probe, indices[probe]), , drop=FALSE]
        ## middle:
        if (nrow(logRASPCF)>1) {
          for (probe in 2:nrow(logRASPCF)) {
            logRC <- rbind(logRC, logRASPCF[rep(probe, indices[probe]-indices[probe-1]), ])
          }
        }
        ## end:
        if (nrow(lr)>indices[probe]) {
          logRC <- rbind(logRC, logRASPCF[rep(probe, nrow(lr)-indices[probe]), ])
        }

        ## Unlike in the single sample implementation, logR levels do not need to be adjusted as breakpoints haven't changed
        ## Similarly, germline homozygous stretches do not need to be resegmented as they have been included in the segmentation from the start


        ## fill in  NAs in logR data (otherwise they cause problems later on)
        ## some NA probes are filled in with zero, replace those too.
        ## In this implementation we do this for each chrke separately
        for (samplei in 1:ncol(logRC)) {
          nakes = c(which(is.na(logRC[, samplei])), which(logRC[, samplei]==0))
          nonnakes = which(!is.na(logRC[, samplei])&!(logRC[, samplei]==0))
          if (length(nakes)>0) {
            for (nake in 1:length(nakes)) {
              pna = nakes[nake]
              if (length(nonnakes)>0) {
                closestnonna = nonnakes[which.min(abs(nonnakes-pna))]
                logRC[pna, samplei] = logRC[closestnonna, samplei]
              } else {
                logRC[pna, samplei] = 0
              }
            }
          }
        }

        ## combine result with that of previous chrkes
        logRPCFed <- rbind(logRPCFed, logRC)

      } ## end of SEGMENTATION
    } ## end for chrke

    rownames(bafPCFed) = bafnames
    ## remove logRonly Sites from segmented baf data (these sites cause problems with ascat otherwise)

    rownames(logRPCFed) = rownames(ASCATobj$Tumor_LogR)

    #adapt levels again
    for (samplei in 1:length(ASCATobj$samples)) {
      seg = rle(logRPCFed[, samplei])$lengths
      logRPCFed_s = numeric(0)
      startprobe = 1
      endprobe = 0
      prevlevel = 0
      for (i in 1:length(seg)) {
        endprobe = endprobe+seg[i]
        level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe, samplei], na.rm=TRUE)
        #making sure no NA's get filled in...
        if (is.nan(level)) {
          level=prevlevel
        } else {
          prevlevel=level
        }
        logRPCFed_s = c(logRPCFed_s, rep(level, seg[i]))
        startprobe = startprobe + seg[i]
      }
      logRPCFed[, samplei] <- logRPCFed_s
    }

    # if less than 800 segments: this segmentlength is ok, otherwise, rerun with higher segmentlength
    if (all(apply(logRPCFed, 2, function(x) length(unique(x)))<800)) {
      break
    }
  } ## end for segmentlength in segmentlengths

  # add column names and write results to files
  colnames(logRPCFed) = colnames(ASCATobj$Tumor_LogR)
  colnames(bafPCFed) = colnames(ASCATobj$Tumor_LogR)
  Tumor_BAF_segmented <- list()

  for (sample in 1:length(ASCATobj$samples)) {
    logrfilename = paste(out.dir, "/", ASCATobj$samples[sample], ".LogR.PCFed.txt", sep="")
    baffilename = paste(out.dir, "/", ASCATobj$samples[sample], ".BAF.PCFed.txt", sep="")

    ## remove NAs from BAF data
    bafPCFed_sample <- bafPCFed[!is.na(bafPCFed[, sample]), sample, drop=FALSE]
    Tumor_BAF_segmented[[sample]] <- 1-bafPCFed_sample

    if (!is.na(out.dir)) write.table(logRPCFed[, sample], logrfilename, sep="\t", col.names=FALSE, quote=FALSE)
    if (!is.na(out.dir)) write.table(bafPCFed_sample, baffilename, sep="\t", col.names=FALSE, quote=FALSE)
  }

  ASCATobj$Tumor_LogR_segmented=logRPCFed
  ASCATobj$Tumor_BAF_segmented=Tumor_BAF_segmented
  ASCATobj$failedarrays=ascat.gg$failedarrays
  return(ASCATobj)

}


################################################################################
## helper functions

runFastASMultiPCF <- function(x, w, gamma, yest, subsize=5000, wsample, nProbesPerSeg) {
  ## Runs exact algorithm on overlapping windows and save breakpoints.
  ## Then runs compact algorithm globally on all collected breakpoints and
  ## Takes samples in columns
  ## Adapted from runMultiPcfSubset (copynumber package)
  antGen <- nrow(x)
  SUBSIZE <- subsize   #length of subsets
  mark <- rep(0, antGen)
  start0=1
  while (start0 + SUBSIZE < antGen) {
    ## segment part of chromosome and save breakpoints
    res <- ASmultiPCFcompact(nr=w[, start0 + (0:SUBSIZE)], wSum=t(x[start0 + (0:SUBSIZE), ]), gamma=gamma, wsample=wsample, nProbesPerSeg = rep(1, SUBSIZE+1))
    mark[start0+res$sta-1] <- 1
    ## find start and end of new chromosome part
    start0 <- start0 + 4*SUBSIZE/5
  }
  ## add all remaining data points to break points
  mark[start0:antGen] <- 1

  ## compress data based on break points
  compX <- compactASMulti(y=t(x), mark = as.logical(mark), w = w)
  ## run exact algorithm on compressed data
  compPotts <- ASmultiPCFcompact(nr=compX$Weight, wSum=compX$wSum, gamma = gamma, wsample=wsample, nProbesPerSeg=compX$nProbesPerSeg)

  if (yest) {
    potts <- expandMulti(nrow(x), ncol(x), compPotts$Lengde, compPotts$mean)
    return(list(yhat = potts, length = compPotts$Lengde, start0 = compPotts$sta,
                mean = compPotts$mean, nIntervals = compPotts$nIntervals))
  } else {
    return(list(length = compPotts$Lengde, start0 = compPotts$sta,
                mean = compPotts$mean, nIntervals = compPotts$nIntervals))
  }

}


# function that accumulates numbers of observations, weight of observations and
# weighted sums between potential breakpoints
compactASMulti <- function(y, mark, w) {
  antGen <- ncol(y)
  antSample <- nrow(y)
  antMark <- sum(mark)
  ant <- rep(0, antMark)
  sum <- weight <- matrix(0, nrow=antSample, ncol=antMark)
  nProbesPerSeg <- rep(0, antMark)
  pos <- 1
  oldPos <- 0
  count <- 1
  delSum <- delWeight <- rep(0, antSample)
  while (pos <= antGen) {
    delSum <- delWeight <- rep(0, antSample)
    while (mark[pos] < 1) {
      delSum <- delSum + y[, pos]*w[, pos] ## weighted sum
      delWeight <- delWeight + w[, pos]
      pos <- pos+1
    }
    ant[count] <- pos-oldPos
    sum[, count] <- delSum+y[, pos]*w[, pos] ## weighted sum
    weight[, count] <- delWeight+w[, pos]
    if (count==1) {
      nProbesPerSeg[count] <- pos
    } else {
      nProbesPerSeg[count] <- pos - sum(nProbesPerSeg[1:(count-1)])
    }
    oldPos <- pos
    pos <- pos+1
    count <- count+1
  }
  list(wSum=sum, Weight=weight, nProbesPerSeg=nProbesPerSeg)
}


## main calculations for fast allele-specific multipcf-versions
ASmultiPCFcompact <- function(nr, wSum, gamma, wsample, nProbesPerSeg) {
  ## nr, wSum : numbers and weighted(!) sums for one analysis unit,
  ## typically one chromosomal arm. Samples assumed to be in rows.
  ## gamma: penalty for discontinuities
  N <- ncol(nr)
  nSamples <- nrow(wSum)
  ## initialisations
  yhat <- rep(0, N*nSamples)
  dim(yhat) <- c(nSamples, N)
  bestCost <- rep(0, N)
  bestSplit <- rep(0, N+1)
  bestAver <- rep(0, N*nSamples)
  dim(bestAver) <- c(nSamples, N)
  Sum <- rep(0, N*nSamples)
  dim(Sum) <- c(nSamples, N)
  Nevner <- rep(0, N*nSamples)
  dim(Nevner) <- c(nSamples, N)
  eachCost <- rep(0, N*nSamples)
  dim(eachCost) <- c(nSamples, N)
  Cost <- rep(0, N)
  ## Filling of first elements
  Sum[, 1]<-wSum[, 1]
  Nevner[, 1]<-nr[, 1]
  bestSplit[1]<-0
  bestAver[, 1] <- wSum[, 1]/nr[, 1]
  helper <- rep(1, nSamples)
  bestCost[1]<-helper %*% (-Sum[, 1]*bestAver[, 1])
  lengde <- rep(0, N)

  ## Solving for gradually longer arrays. Sum accumulates
  ## error values for righthand plateau downward from n
  ## this error added to gamma and the stored cost in bestCost
  ## give the total cost. Cost stores the total cost for breaks
  ## at any position below n, and which.min finds the position
  ## with lowest cost (best split). Aver is the average of the
  ## righthand plateau.
  for (n in 2:N) {
    Sum[, 1:n] <- Sum[, 1:n]+wSum[, n] ## wSum is already weighted!
    Nevner[, 1:n] <- Nevner[, 1:n]+nr[, n]
    if (any(Nevner[, 1:n]==0)) {
      Nevnertmp <- Nevner[, 1:n]
      Nevnertmp[Nevnertmp==0] <- 1
      eachCost[, 1:n] <- -(Sum[, 1:n]^2)/Nevnertmp
    } else {
      eachCost[, 1:n] <- -(Sum[, 1:n]^2)/Nevner[, 1:n]
    }
    Cost[1:n] <- helper %*% eachCost[, 1:n]
    Cost[2:n] <- Cost[2:n]+bestCost[1:(n-1)]+gamma
    Pos <- which.min(Cost[1:n])
    bestCost[n] <- Cost[Pos]
    bestAver[, n] <- Sum[, Pos]/Nevner[, Pos]
    bestSplit[n] <- Pos-1
  }

  ## The final solution is found iteratively from the sequence
  ## of split positions stored in bestSplit and the averages
  ## for each plateau stored in bestAver

  n <- N
  antInt <- 0
  while (n > 0) {
    yhat[, (bestSplit[n]+1):n] <- bestAver[, n]
    antInt <- antInt+1
    lengde[antInt] <- sum(nProbesPerSeg[(bestSplit[n]+1):n])
    n <- bestSplit[n]
  }
  lengdeRev <- lengde[antInt:1]
  init <- rep(0, antInt)
  init[1]<-1
  if (antInt>=2) {
    for (k in 2:antInt) {
      init[k]<-init[k-1]+lengdeRev[k-1]
    }
  }

  n <- N
  verdi <- rep(0, antInt*nSamples)
  dim(verdi) <- c(nSamples, antInt)
  bestSplit[n+1] <- n
  antall <- antInt
  while (n > 0) {
    verdi[, antall] <- bestAver[, n]
    n <- bestSplit[n]
    antall <- antall-1
  }

  list(Lengde = lengdeRev, sta = init, mean = verdi, nIntervals=antInt)
}

## expand compact solution
expandMulti <- function(nProbes, nSamples, lengthInt, mean) {
  ##input: nr of probes, length of intervals,
  ## value in intervals; returns the expansion

  Potts <- rep(0, nProbes*nSamples)
  dim(Potts) <- c(nSamples, nProbes)
  lengthCompArr <- length(lengthInt)
  k <- 1
  for (i in 1:lengthCompArr) {
    for (j in 1:lengthInt[i]) {
      Potts[, k] <- mean[, i]
      k <- k+1
    }
  }
  return(Potts)
}

################################################################################
## matrix version of basic ASCAT functions

## mirror BAF values around 0.5 to reduce regions with allelic imbalance to one band
## (needed for winsorization and segmentation)
mirrorBafMatrix <- function(baf) {
  ## samples in columns
  apply(baf, 1:2, function(x) {
    ifelse(x>0.5, x, 1-x)
  })
}

## winsorize a matrix of samples on a per sample base while taking care of NAs
madWinsMatrixWithNA <- function(x, tau, k) {
  ## samples in columns
  if (nrow(x)==1) {
    xwins <- x
    xwins[!is.na(x)] = madWins(x[!is.na(x)], tau, k)$ywin
    return(xwins)
  } else {
    apply(x, 2, function(x) {
      xwins <- vector(mode="numeric", length=length(x))
      xwins[is.na(x)] = NA
      xwins[!is.na(x)] = madWins(x[!is.na(x)], tau, k)$ywin
      return(xwins)
    })
  }
}


#Get mad SD (based on KL code)
getMadwithNA <- function(x, k=25) {

  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0&!is.na(x)]

  #Calculate runMedian
  runMedian <- medianFilter(x, k)

  dif <- x-runMedian
  SD <- mad(dif)

  return(SD)
}
