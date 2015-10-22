# function to read in SNP array data
# input: filenames of tumor LogR, tumor BAF, germline LogR and germline BAF data
# germline data files can be NULL - in that case these are not read in
# output: a data structure containing:
# 1. Tumor_LogR data matrix
# 2. Tumor_BAF data matrix
# 3. Tumor_LogR_segmented: placeholder, NULL
# 4. Tumor_BAF_segmented: placeholder, NULL
# 5. Germline_LogR data matrix
# 6. Germline_BAF data matrix
# 7. SNPpos: position of all SNPs
# 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13
# 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)
ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X")) {

  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")

  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA

  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")

    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }

  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]

  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  }
 
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))  
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }

  # write final SNP positions to file
  #write.table(SNPpos,"SNPpos.txt", sep="\t")

  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chuncks...")
  chr = split_genome(SNPpos)

  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
              samples = colnames(Tumor_LogR)))
}




# plots SNP array data
# input: an ASCAT object (e.g. from ASCAT.loaddata and plots the SNP array data
# tumorfiles: start of filename for tumor data plots (no plotting if NULL)
# germlinefiles: start of filename for germline data plots (no plotting if NULL)
ascat.plotRawData = function(ASCATobj, tumorfiles = "Tumor", germlinefiles = "Germline") {
  attach(ASCATobj)
  if(!is.null(tumorfiles)) {
    print.noquote("Plotting tumor data")
    for (i in 1:dim(Tumor_LogR)[2]) {
      png(filename = paste(tumorfiles,samples[i],".png",sep=""), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(Tumor_LogR)[1]>20000,".",20))
      plot(c(1,dim(Tumor_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
      points(Tumor_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
      points(Tumor_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
  if(!is.null(germlinefiles) && !is.null(Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(Germline_LogR)[2]) {
      png(filename = paste(germlinefiles,samples[i],".png",sep=""), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(Tumor_LogR)[1]>20000,".",20))
      plot(c(1,dim(Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(Germline_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
  detach(ASCATobj)
}



# run ASPCF segmentation
# input: (i) an ASCAT object from e.g. ASCAT.loaddata; 
# (ii) selectsamples: a vector containing the sample(number)s to PCF otherwise
# (iii) germline genotypes (NULL if germline data is available)
# output: a data structure containing:
# 1. Tumor_LogR data matrix
# 2. Tumor_BAF data matrix
# 3. Tumor_LogR_segmented: matrix of LogR segmented values
# 4. Tumor_BAF_segmented: list of BAF segmented values; each element in the list is a matrix containing the segmented values for one sample (only for probes that are germline homozygous)
# 5. Germline_LogR data matrix
# 6. Germline_BAF data matrix
# 7. SNPpos: position of all SNPs
# 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13
# 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)
# note: this function can be easily parallelized by controlling the selectsamples parameter
# it saves the results in LogR_PCFed[sample]_[segment].txt and BAF_PCFed[sample]_[segment].txt
# if these files exist, the results are read from the files
# hence, after parallelization, copy all the files into the current directory, and run this function to read in the results
ascat.aspcf = function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL) {
  attach(ASCATobj)
 
  #first, set germline genotypes
  gg = NULL
  if(!is.null(ascat.gg)) {
    gg = ascat.gg
  }
  else {
    gg = Germline_BAF < 0.3 | Germline_BAF > 0.7
  }
  # calculate germline homozygous stretches for later resegmentation
  ghs = predictGermlineHomozygousStretches(chr, gg)

  segmentlengths = c(25,50,100,200,400,800)

  Tumor_LogR_segmented = matrix(nrow = dim(Tumor_LogR)[1], ncol = dim(Tumor_LogR)[2])
  rownames(Tumor_LogR_segmented) = rownames(Tumor_LogR)
  colnames(Tumor_LogR_segmented) = colnames(Tumor_LogR)
  Tumor_BAF_segmented = list();
  source("aspcf.R")
  for (sample in selectsamples) {
    print.noquote(paste("Sample ", samples[sample], " (",sample,"/",length(samples),")",sep=""))
    logrfilename = paste("LogR_PCFed_",samples[sample],".txt",sep="")
    baffilename = paste("BAF_PCFed_",samples[sample],".txt",sep="")
    logRPCFed = numeric(0)
    bafPCFed = numeric(0)
    if(length(dir(pattern=logrfilename))==0 || length(dir(pattern=baffilename))==0) {
      for (segmentlength in segmentlengths) {
        logRPCFed = numeric(0)
        bafPCFed = numeric(0)
        tbsam = Tumor_BAF[,sample]
        names(tbsam) = rownames(Tumor_BAF)
        homosam = gg[,sample]
        for (chrke in 1:length(chr)) {
          lr = Tumor_LogR[chr[[chrke]],sample]
          #winsorize to remove outliers
          #this has a problem with NAs
          lrwins = vector(mode="numeric",length=length(lr))
          lrwins[is.na(lr)] = NA
	  lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)],2.5,25)$ywin
          baf = tbsam[chr[[chrke]]]
          homo = homosam[chr[[chrke]]]
          Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
          bafsel = baf[Select_het]
          indices = which(Select_het)
          logRaveraged = NULL;
          if(length(indices)!=0) {
            averageIndices = c(1,(indices[1:(length(indices)-1)]+indices[2:length(indices)])/2,length(lr)+0.01)
            startindices = ceiling(averageIndices[1:(length(averageIndices)-1)])
            endindices = floor(averageIndices[2:length(averageIndices)]-0.01)
            if(length(indices)==1) {
              startindices = 1
              endindices = length(lr)
            }
            nrIndices = endindices - startindices + 1
            logRaveraged = vector(mode="numeric",length=length(indices))
            for(i in 1:length(indices)) {
              if(is.na(endindices[i])) {
                endindices[i]=startindices[i]
              }
              logRaveraged[i]=mean(lrwins[startindices[i]:endindices[i]], na.rm=T)
            }
          }
          # if there are no probes in the segment (after germline homozygous removal), don't do anything, except add a LogR segment
          if(length(logRaveraged)>0) {
            logRASPCF = NULL
            bafASPCF = NULL
            if(length(logRaveraged)<6) {
              logRASPCF = rep(mean(logRaveraged),length(logRaveraged))
              bafASPCF = rep(mean(bafsel),length(logRaveraged))
            }
            else {
              PCFed = fastAspcf(logRaveraged,bafsel,3,segmentlength)
              logRASPCF = PCFed$yhat1
              bafASPCF = PCFed$yhat2
            }
            names(bafASPCF)=names(indices)
            logRc = numeric(0)
            for(probe in 1:length(logRASPCF)) {
              if(probe == 1) {
                logRc = rep(logRASPCF[probe],indices[probe])
              }
              # if probe is 1, set the beginning, and let the loop go:
              if(probe == length(logRASPCF)) {
                logRc = c(logRc,rep(logRASPCF[probe],length(lr)-indices[probe]))
              }
              else if(logRASPCF[probe]==logRASPCF[probe+1]) {
                logRc = c(logRc, rep(logRASPCF[probe],indices[probe+1]-indices[probe]))
              }
              else {
                #find best breakpoint
                d = numeric(0)
                totall = indices[probe+1]-indices[probe]
                for (bp in 0:(totall-1)) {
                  dis = sum(abs(lr[(1:bp)+indices[probe]]-logRASPCF[probe]), na.rm=T)
                  if(bp!=totall) {
                    dis = sum(dis, sum(abs(lr[((bp+1):totall)+indices[probe]]-logRASPCF[probe+1]), na.rm=T), na.rm=T)
                  }
                  d = c(d,dis)  
                }
                breakpoint = which.min(d)-1
                logRc = c(logRc,rep(logRASPCF[probe],breakpoint),rep(logRASPCF[probe+1],totall-breakpoint))
              }
            }
            #2nd step: adapt levels!
            logRd = numeric(0)
            seg = make_seg_lr(logRc)
            startprobe = 1
            endprobe = 0
            for (i in 1:length(seg)) {
              endprobe = endprobe+seg[i]
              level = mean(lr[startprobe:endprobe], na.rm=T)
              logRd = c(logRd, rep(level,seg[i]))
              startprobe = startprobe + seg[i]
            }
            logRPCFed = c(logRPCFed,logRd)
            bafPCFed = c(bafPCFed,bafASPCF)
          }
          # add a LogR segment
          else {
            level = mean(lr,na.rm=T)
            reps = length(lr)
            logRPCFed = c(logRPCFed,rep(level,reps)) 
          }
          # correct wrong segments in germline homozygous stretches:
          homsegs = ghs[[sample]][ghs[[sample]][,1]==chrke,]  
          startchr = min(chr[[chrke]])
          endchr = max(chr[[chrke]])
          # to solve an annoying error when homsegs has length 1:
          if(length(homsegs)==3) {
            homsegs=t(as.matrix(homsegs))
          }
          if(!is.null(homsegs)&&!is.na(homsegs)&&dim(homsegs)[1]!=0) {
            for (i in 1:dim(homsegs)[1]) {
              # note that only the germline homozygous segment is resegmented, plus a bit extra (but this is NOT replaced)
              startpos = max(homsegs[i,2],startchr)
              endpos = min(homsegs[i,3],endchr)
              # PCF over a larger fragment
              startpos2 = max(homsegs[i,2]-100,startchr)
              endpos2 = min(homsegs[i,3]+100,endchr)
              # take into account a little extra (difference between startpos2 and startpos3 is not changed)
              startpos3 = max(homsegs[i,2]-5,startchr)
              endpos3 = min(homsegs[i,3]+5,endchr)
              # note that the parameters are arbitrary, but <100 seems to work on the ERBB2 example!
              # segmentlength is lower here, as in the full data, noise on LogR is higher!
              # do this on Winsorized data too!
              towins = Tumor_LogR[startpos2:endpos2,sample]
              winsed = madWins(towins[!is.na(towins)],2.5,25)$ywin
              pcfed = vector(mode="numeric",length=length(towins))
              pcfed[!is.na(towins)] = exactPcf(winsed,6,floor(segmentlength/4))
              pcfed2 = pcfed[(startpos3-startpos2+1):(endpos3-startpos2+1)]
              dif = abs(pcfed2-logRPCFed[startpos3:endpos3])
              #0.3 is hardcoded here, in order not to have too many segments!
              #only replace if enough probes differ (in order not to get singular probes with different breakpoints)
              if(!is.na(dif)&&sum(dif>0.3)>5) {
                #replace a bit more to make sure no 'lone' probes are left (startpos3 instead of startpos)
                logRPCFed[startpos3:endpos3]=ifelse(dif>0.3,pcfed2,logRPCFed[startpos3:endpos3])
              }
            }
          }
        }
        #fill in NAs (otherwise they cause problems):
        #some NA probes are filled in with zero, replace those too:
        nakes = c(which(is.na(logRPCFed)),which(logRPCFed==0))
        nonnakes = which(!is.na(logRPCFed)&!(logRPCFed==0))
        if(length(nakes)>0) {
          for (nake in 1:length(nakes)) {
            pna = nakes[nake]
            closestnonna = which.min(abs(nonnakes-pna))
            logRPCFed[pna] = logRPCFed[closestnonna]
          }
        }
        #adapt levels again
        seg = make_seg_lr(logRPCFed)
        logRPCFed = numeric(0)
        startprobe = 1
        endprobe = 0
        prevlevel = 0
        for (i in 1:length(seg)) {
          endprobe = endprobe+seg[i]
          level = mean(Tumor_LogR[startprobe:endprobe,sample], na.rm=T)
          #making sure no NA's get filled in...
          if(is.nan(level)) {
            level=prevlevel
          }
          else {
            prevlevel=level
          }
          logRPCFed = c(logRPCFed, rep(level,seg[i]))
          startprobe = startprobe + seg[i]
        }
        #put in names and write results to files
        names(logRPCFed) = rownames(Tumor_LogR)

        # if less than 800 segments: this segmentlength is ok, otherwise, rerun with higher segmentlength
        if(length(unique(logRPCFed))<800) {
          break
        }
      }

      write.table(logRPCFed,logrfilename,sep="\t",col.names=F)
      write.table(bafPCFed,baffilename,sep="\t",col.names=F)
      bafPCFed = as.matrix(bafPCFed)
    }
    else {
      logRPCFed = read.table(logrfilename, sep="\t", header=F, row.names=1)[,1]
      bafPCFed = read.table(baffilename, sep="\t", header=F, row.names=1)
    }
    Tumor_LogR_segmented[,sample] = logRPCFed
    Tumor_BAF_segmented[[sample]] = 1-bafPCFed
  }
  
  
  ASCATobj = list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
                  Tumor_LogR_segmented = Tumor_LogR_segmented, 
                  Tumor_BAF_segmented = Tumor_BAF_segmented, 
                  Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
                  SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
                  samples = colnames(Tumor_LogR))
  detach(ASCATobj)
  return(ASCATobj)
}



# plots SNP array data
# input: an ASCAT object (e.g. from ASCAT.ASPCF) and plots the SNP array data before and after segmentation
# filenames: start of the names of the output files
ascat.plotSegmentedData = function(ASCATobj, filenames = "ASPCF") {
  attach(ASCATobj)
  for (arraynr in 1:dim(Tumor_LogR)[2]) {
    Select_nonNAs = rownames(Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(Tumor_LogR)[1]
    names(AllIDs) = rownames(Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = paste(filenames,samples[arraynr],".png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2)
    r = Tumor_LogR_segmented[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr]
    beta = Tumor_BAF_segmented[[arraynr]][,]
    plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
    points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
    points(r,col="green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = intersect(ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
    points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
    points(beta, col = "green")
    points(1-beta, col = "green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = intersect(ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
  detach(ASCATobj)
}


# the ASCAT main function, calculating the allele-specific copy numbers
# input: (i) an ASCAT object from ASCAT.ASPCF.
# (ii) gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
# (iii) sunrisefiles: filenames of 'sunrise plots' (plots with distance matrices)
# (iv) profilefiles: filenames of ASCAT profiles (plots with final allele-specific copy number)
# (v) rawprofilefiles: filenames of raw copy number profiles
# (vi) allow100percent: setting thris to TRUE allows solutions with 100 % aberrant cell fraction, e.g. cell lines, non-tumor tissue (with flat profiles) - not advisably unless you have a clear reason to do so
# output: an ASCAT output object, containing:
# 1. nA: copy number of the A allele
# 2. nB: copy number of the B allele
# 3. aberrantcellfraction: the aberrant cell fraction of all arrays
# 4. ploidy: the ploidy of all arrays
# 5. failedarrays: arrays on which ASCAT analysis failed
# 6. segments: a list containing the copy number segments of each sample (with NA for failed arrays)
# note: for copy number only probes, nA contains the copy number value and nB = 0.
ascat.runAscat = function(ASCATobj, gamma = 0.55, sunrisefiles = "sunrise", profilefiles = "ASCATprofile", rawprofilefiles = "rawprofile", allow100percent = F) {
  attach(ASCATobj)
  goodarrays=NULL
  res = vector("list",dim(Tumor_LogR)[2])
  for (arraynr in 1:dim(Tumor_LogR)[2]) {
    print.noquote(paste("Sample ", samples[arraynr], " (",arraynr,"/",length(samples),")",sep=""))
    lrr=Tumor_LogR[,arraynr]
    names(lrr)=rownames(Tumor_LogR)
    baf=Tumor_BAF[,arraynr]
    names(baf)=rownames(Tumor_BAF)
    lrrsegm = Tumor_LogR_segmented[,arraynr]
    names(lrrsegm) = rownames(Tumor_LogR_segmented)
    bafsegm = Tumor_BAF_segmented[[arraynr]][,]
    names(bafsegm) = rownames(Tumor_BAF_segmented[[arraynr]])
    res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ch,paste(sunrisefiles,samples[arraynr],".png",sep=""),paste(profilefiles,samples[arraynr],".png",sep=""),paste(rawprofilefiles,samples[arraynr],".png",sep=""),gamma,allow100percent)
    if(!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays)+1] = arraynr
    }
  }

  if(length(goodarrays)>0) {
    n1 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
    n2 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) = rownames(Tumor_LogR)
    rownames(n2) = rownames(Tumor_LogR)
    colnames(n1) = colnames(Tumor_LogR)[goodarrays]
    colnames(n2) = colnames(Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[,i] = res[[goodarrays[i]]]$nA
      n2[,i] = res[[goodarrays[i]]]$nB
    }

    tp = vector(length=length(goodarrays))
    psi = vector(length=length(goodarrays))
    ploidy = vector(length=length(goodarrays))
    for (i in 1:length(goodarrays)) {
      tp[i] = res[[goodarrays[i]]]$rho
      psi[i] = res[[goodarrays[i]]]$psi
      ploidy[i] = mean(res[[goodarrays[i]]]$nA+res[[goodarrays[i]]]$nB,na.rm=T)
    }
    fa = colnames(Tumor_LogR)[-goodarrays]
  }
  else {
    n1 = NULL
    n2 = NULL
    tp = NULL
    ploidy = NULL
    psi = NULL
    fa = colnames(Tumor_LogR)
  }

  seg = list()
  for (i in 1:length(samples)) {
    seg[[i]] = res[[i]]$seg
  }

  detach(ASCATobj)

  return(list(nA = n1, nB = n2, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, failedarrays = fa, segments = seg))
}





# helping function to read segments:
make_seg_lr = function(r) {
  pcf_segments = numeric(0);
  index = 0;
  previousr = 1E10;
  for (i in 1:length(r)) {
    if (r[i] != previousr) {
      index=index+1;
      count=1;
    }
    else {
      count = count + 1;
    }
    pcf_segments[index] = count;
    previousr = r[i];
  }
  return(pcf_segments);
}



# helper function to split the genome into parts
split_genome = function(SNPpos) {

  # look for gaps of more than 1Mb and chromosome borders
  holesOver1Mb = which(diff(SNPpos[,2])>=1000000)+1
  chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1),1]!=SNPpos[2:(dim(SNPpos)[1]),1])+1

  holes = unique(sort(c(holesOver1Mb,chrBorders)))

  # find which segments are too small
  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  while (1 %in% joincandidates) {
    holes=holes[-1]
    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
  while ((length(holes)+1) %in% joincandidates) {
    holes=holes[-length(holes)]
    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
 
  while(length(joincandidates)!=0) {
    # the while loop is because after joining, segments may still be too small..

    startseg = c(1,holes)
    endseg = c(holes-1,dim(SNPpos)[1])

    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    previoussamechr = SNPpos[endseg[joincandidates-1],1]==SNPpos[startseg[joincandidates],1] 
    nextsamechr = SNPpos[endseg[joincandidates],1]==SNPpos[startseg[joincandidates+1],1]

    distanceprevious = SNPpos[startseg[joincandidates],2]-SNPpos[endseg[joincandidates-1],2]
    distancenext = SNPpos[startseg[joincandidates+1],2]-SNPpos[endseg[joincandidates],2]

    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    joins = ifelse(previoussamechr&nextsamechr, 
                   ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
                   ifelse(nextsamechr, joincandidates, joincandidates-1))

    holes=holes[-joins]

    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.

  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPpos)[1])

  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }
  
  return(chr)
}



# helper function to predict germline homozygous segments for later resegmentation
predictGermlineHomozygousStretches = function(chr, hom) {

  # contains the result: a list of vectors of probe numbers in homozygous stretches for each sample
  HomoStretches = list()

  for (sam in 1:dim(hom)[2]) {
    homsam = hom[,sam]

    perchom = sum(homsam,na.rm=T)/sum(!is.na(homsam))

    # NOTE THAT A P-VALUE THRESHOLD OF 0.001 IS HARDCODED HERE
    homthres = ceiling(log(0.001,perchom))

    allhprobes = NULL
    for (chrke in 1:length(chr)) {
      hschr = homsam[chr[[chrke]]]

      hprobes = vector(length=0)
      for(probe in 1:length(hschr)) {
        if(!is.na(hschr[probe])) {
          if(hschr[probe]) {
            hprobes = c(hprobes,probe)
          }
          else {
            if(length(hprobes)>=homthres) {
              allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
            }
            hprobes = vector(length=0)
          }
        }
      }
      # if the last probe is homozygous, this is not yet accounted for
      if(!is.na(hschr[probe]) & hschr[probe]) {
        if(length(hprobes)>=homthres) {
          allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
        }
      }
   
    }

    HomoStretches[[sam]]=allhprobes
  
  }

  return(HomoStretches)
}





# function to make segments of constant LRR and BAF 
# this function is more general and does not depend on specifically ASPCF output
# it can also handle segmention performed on LRR and BAF separately
make_segments = function(r,b) {
  m = matrix(ncol = 2, nrow = length(b))
  m[,1] = r
  m[,2] = b
  m = as.matrix(na.omit(m))
  pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) = c("r","b","length");
  index = 0;
  previousb = -1;
  previousr = 1E10;
  for (i in 1:dim(m)[1]) {
    if (m[i,2] != previousb || m[i,1] != previousr) {
      index=index+1;
      count=1;
      pcf_segments[index, "r"] = m[i,1];
      pcf_segments[index, "b"] = m[i,2];
    }
    else {
      count = count + 1;
    }
    pcf_segments[index, "length"] = count;
    previousb = m[i,2];
    previousr = m[i,1];
  }
  pcf_segments = as.matrix(na.omit(pcf_segments))[,]
  return(pcf_segments);
}



# function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
# input: segmented LRR and BAF and the value for gamma
create_distance_matrix = function(segments, gamma) {
  s = segments
  psi_pos = seq(1,6,0.05) 
  rho_pos = seq(0.1,1.05,0.01)
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      # choose the minor allele
      nMinor = NULL
      if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
        nMinor = nA
      }
      else {
        nMinor = nB
      }
      d[i,j] = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T) 
    }
  }
  return(d)
}



# the ASCAT main function
# lrr: (unsegmented) log R, in genomic sequence (all probes), with probe IDs
# baf: (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
# lrrsegmented: log R, segmented, in genomic sequence (all probes), with probe IDs
# bafsegmented: B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
# gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
# chromosomes: a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
# distancepng: if NA: distance is plotted, if filename is given, the plot is written to a .png file
# copynumberprofilespng: if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file
# nonroundedprofilepng: if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, gamma = 0.55, allow100percent) {
  ch = chromosomes
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  library(RColorBrewer)

  s = make_segments(r,b)
  d = create_distance_matrix(s, gamma)

  if (!is.na(distancepng)) {
    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
  }

  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)

  hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")

  axis(1, at = seq(0, 1, by = 1/5), label = seq(1, 6, by = 1))
  axis(2, at = seq(0, 1/1.05, by = 1/3/1.05), label = seq(0.1, 1, by = 3/10))

  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)

  nropt = 0
  localmin = NULL
  optima = list()
  for (i in 4:(dim(d)[1]-3)) {
    for (j in 4:(dim(d)[2]-3)) {
      m = d[i,j]
      seld = d[(i-3):(i+3),(j-3):(j+3)]
      seld[4,4] = max(seld)
      if(min(seld) > m) {
        psi = as.numeric(rownames(d)[i])
        rho = as.numeric(colnames(d)[j])
        nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        
        # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
        ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
        percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
        perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
        # the next can happen if BAF is a flat line at 0.5 
        if (is.na(perczeroAbb)) {
          perczeroAbb = 0
        }

        goodnessOfFit = (1-m/TheoretMaxdist) * 100

        if (ploidy > 1.5 & ploidy < 5.5 & rho >= 0.2 & goodnessOfFit > 75 & (percentzero > 0.01 | perczeroAbb > 0.1)) {
          nropt = nropt + 1
          optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
          localmin[nropt] = m
        }
      }     
    }
  }
  # if solutions with 100 % aberrant cell fraction should be allowed:
  # if there are no solutions, drop the conditions on regions with copy number zero, and include the borders (rho = 1) as well
  # this way, if there is another solution, this is still preferred, but these solutions aren't standardly eliminated
  if (allow100percent & nropt == 0) {
    #first, include borders
    cold = which(as.numeric(colnames(d))>1)
    d[,cold]=1E20
    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i,j]
        seld = d[(i-3):(i+3),(j-3):(j+3)]
        seld[4,4] = max(seld)
        if(min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])
          nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
          nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        
          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
          percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
          perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
          # the next can happen if BAF is a flat line at 0.5 
          if (is.na(perczeroAbb)) {
            perczeroAbb = 0
          }

          goodnessOfFit = (1-m/TheoretMaxdist) * 100


          if (ploidy > 1.5 & ploidy < 5.5 & rho >= 0.2 & goodnessOfFit > 75) {
            nropt = nropt + 1
            optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }
  }

  if (nropt>0) {
    optlim = sort(localmin)[1]
    for (i in 1:length(optima)) {
      if(optima[[i]][1] == optlim) {  
        psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
        rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
        if(rho_opt1 > 1) {
          rho_opt1 = 1
        }
        ploidy_opt1 = optima[[i]][4]
        goodnessOfFit_opt1 = optima[[i]][5]
        points((psi_opt1-1)/5,(rho_opt1-0.1)/0.95,col="green",pch="X", cex = 2)
      }
    }
  }

  if (!is.na(distancepng)) {
    dev.off()
  }


  if(nropt>0) {

    if (!is.na(nonroundedprofilepng)) {
      png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200)
    } 
    else {      
      windows(10,5)
    }

    par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)

    rho = rho_opt1
    psi = psi_opt1
    
    nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
    plot(c(1,length(nAfull)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nBfull,col="blue",pch = "|")
    points(nAfull+nBfull,col="purple",pch = "|")
# don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,5,chrs[i], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    
    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }


    rho = rho_opt1
    psi = psi_opt1

#note that this approach may be a bit of a problem, as it assumes that every segment has a unique LogR value..
    tlr = unique(lrrsegmented)
    seg = NULL
    for (i in 1:length(tlr)) {
      logR = tlr[i]
      pr = which(lrrsegmented==logR)
      start = min(pr)
      end = max(pr)
      bafke = bafsegmented[intersect(names(lrrsegmented)[pr],names(bafsegmented))][1]
      #if bafke is NA, this means that we are dealing with a germline homozygous stretch with a copy number change within it.
      #in this case, nA and nB are irrelevant, just their sum matters
      if(is.na(bafke)) {
        bafke=0
      }
      nAraw = (rho-1-(bafke-1)*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho
      nBraw = (rho-1+bafke*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho
# correct for negative values:
      if (nAraw+nBraw<0) {
        nAraw = 0
        nBraw = 0
      }
      else if (nAraw<0) {
        nBraw = nAraw+nBraw
        nAraw = 0  
      }
      else if (nBraw<0) {
        nAraw = nAraw+nBraw
        nBraw = 0  
      }
      # when evidence for odd copy number in segments of BAF = 0.5, assume a deviation.. 
      limitround = 0.5
      nA = ifelse(bafke==0.5,
             ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
               round(nAraw)+1,
               ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround, 
                 round(nAraw)-1,
                 round(nAraw))),
             round(nAraw))
      nB = round(nBraw)
      if (is.null(seg)) {
        seg = t(as.matrix(c(start,end,nA,nB)))
      }
      else {
        seg = rbind(seg,c(start,end,nA,nB))
      }
    }
    colnames(seg)=c("start","end","nA","nB")

    # build helping vector
    chrhelp = vector(length=length(lrrsegmented))
    for (chrnr in 1:length(ch)) {
      chrke = ch[[chrnr]]
      chrhelp[chrke] = chrnr
    }    

    # every repeat joins 2 ends. 20 repeats will join about 1 million ends..
    for (rep in 1:20) {
      seg2=seg
      seg = NULL
      skipnext = F
      for(i in 1:dim(seg2)[1]) {
        if(!skipnext) {
          if(i != dim(seg2)[1] && seg2[i,"nA"]==seg2[i+1,"nA"] && seg2[i,"nB"]==seg2[i+1,"nB"] && 
              chrhelp[seg2[i,"end"]]==chrhelp[seg2[i+1,"start"]]) {
            segline = c(seg2[i,"start"],seg2[i+1,"end"],seg2[i,3:4])
            skipnext = T
          }
          else {
            segline = seg2[i,]
          }

          if (is.null(seg)) {
            seg = t(as.matrix(segline))
          }
          else {
            seg = rbind(seg,segline)
          }
        }
        else {
          skipnext = F
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
      nMajor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nA"]
      nMinor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nB"]
    }

    n1all = vector(length = length(lrrsegmented))
    names(n1all) = names(lrrsegmented)
    n2all = vector(length = length(lrrsegmented))
    names(n2all) = names(lrrsegmented)

    # note: any of these can have length 0
    NAprobes = which(is.na(lrr))
    heteroprobes = setdiff(which(names(lrrsegmented)%in%names(bafsegmented)),NAprobes)
    homoprobes = setdiff(setdiff(which(!is.na(baf)),heteroprobes),NAprobes)
    CNprobes = setdiff(which(is.na(baf)),NAprobes)
 
    n1all[NAprobes] = NA
    n2all[NAprobes] = NA
    n1all[CNprobes] = nMajor[CNprobes]+nMinor[CNprobes]
    n2all[CNprobes] = 0
    heteroprobes2 = names(lrrsegmented)[heteroprobes]
    n1all[heteroprobes] = ifelse(baf[heteroprobes2]<=0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n2all[heteroprobes] = ifelse(baf[heteroprobes2]>0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n1all[homoprobes] = ifelse(baf[homoprobes]<=0.5,nMajor[homoprobes]+nMinor[homoprobes],0)
    n2all[homoprobes] = ifelse(baf[homoprobes]>0.5,nMajor[homoprobes]+nMinor[homoprobes],0)


    # plot ASCAT profile
    if (!is.na(copynumberprofilespng)) {
      png(filename = copynumberprofilespng, width = 2000, height = 1000, res = 200)
    } 
    else {      
      windows(10,5)
    }

    par(mar = c(0.5,5,5,0.5), mfrow=c(2,1), cex = 0.4, cex.main=3, cex.axis = 2.5)

    nA2 = n1all[heteroprobes]
    nB2 = n2all[heteroprobes]
    nA = ifelse(nA2>nB2,nA2,nB2)
    nB = ifelse(nA2>nB2,nB2,nA2)
    maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
    plot(c(1,length(nAfull)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nA-0.1,col="red",pch = "|")
    points(nB+0.1,col="green",pch = "|")
    # don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,5,chrs[i], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    
    rBacktransform = gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
    bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
    rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))
    maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
    plot(c(1,length(nAfull)), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(confidence,col="blue",pch = "|")
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,5,chrs[i], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    

    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }

    return(list(rho = rho_opt1, psi = psi_opt1, nA = n1all, nB = n2all, seg = seg))
  }

  else {
    return(list(rho = NA, psi = NA, nA = NA, nB = NA, seg = NA))
  }
  
}


