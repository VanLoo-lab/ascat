#' Function to plot logR and BAF
#' @noRd
plotLogRandBAF=function(Plotdir,SNPpos,LogR,BAF) {
  myLegend=function(ch) {
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = ch[[j]]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      text(tpos,2,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
  }
  stopifnot(identical(rownames(LogR),rownames(BAF)))
  stopifnot(identical(rownames(LogR),rownames(SNPpos)))
  stopifnot(identical(colnames(LogR),colnames(BAF)))
  samples=colnames(LogR)
  chrs=as.vector(unique(SNPpos[,1]))
  last=0
  ch=list()
  SNPorder=vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke=SNPpos[SNPpos[,1]==chrs[i],]
    chrpos=chrke[,2]
    names(chrpos)=rownames(chrke)
    chrpos=sort(chrpos)
    ch[[i]]=(last+1):(last+length(chrpos))
    SNPorder[ch[[i]]]=names(chrpos)
    last=last+length(chrpos)
  }
  SNPpos=SNPpos[SNPorder,]
  LogR=LogR[SNPorder,,drop=F]
  BAF=BAF[SNPorder,,drop=F]
  for (i in 1:dim(LogR)[2]) {
    png(filename=paste(Plotdir,'/',samples[i],".png",sep=""),width=2000,height=1000,res=200)
    par(mar=c(0.5,5,5,0.5),mfrow=c(2,1),cex=0.4,cex.main=3,cex.axis=2,pch=20)
    # plot logR
    plot(c(1,dim(LogR)[1]),c(-2,2),type="n",xaxt="n",main=paste0(samples[i],", LogR"),xlab="",ylab="")
    points(LogR[,i],col="red")
    myLegend(ch)
    # plot mirrored BAF
    plot(c(1,dim(BAF)[1]),c(0,1),type="n",xaxt="n",main=paste0(samples[i],", BAF"),xlab="",ylab="")
    points(ifelse(runif(length(BAF[,i]))<0.5,BAF[,i],1-BAF[,i]),col="red")
    myLegend(ch)
    dev.off()
  }
}

#' Function to compute logR based on counts (normals only)
#' @noRd
computeLogR=function(NORMAL_COUNTS_tot) {
  logR=data.frame(apply(NORMAL_COUNTS_tot,2,function(x) {
    x=log2(x)
    x=x-mean(x)
    return(x)
  }),check.names=F)
  logRref = apply(logR,1,median)
  logR=data.frame(apply(logR,2,function(x) {
    x = x-logRref
    x = x-mean(x)
    return(x)
  }),check.names=F)
  return(logR)
}

#' Generate a cleaned list of SNPs from normal samples.
#'
#' @param Worksheet A table with the following columns: Patient_ID, Normal_ID, Normal_file and Gender.
#' @param Workdir The folder where output should go.
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr").
#' @param minCounts Minimum depth, in normal samples, required for a SNP to be considered.
#' @param X_nonPAR Vector containing genomic coordinates (start & stop) of nonPAR region on X. Default=NULL
#' @param chrom_names A vector containing the names of chromosomes to be considered (optional, default=c(1:22,'X')).
#' @param plotQC A boolean to generate QC reports as PNGs (optional, default=T).
#' @noRd
getLociFromNormals=function(Worksheet, Workdir, alleles.prefix, minCounts, X_nonPAR=NULL, chrom_names=c(1:22,'X'), plotQC=T) {
  stopifnot((is.null(X_nonPAR)) || (length(X_nonPAR)==2 && all(is.numeric(X_nonPAR))))
  # Read all alleleCount files (counts>=0)
  print('      Getting allele counts...')
  NORMAL_COUNTS=foreach(INDEX=1:nrow(Worksheet)) %dopar% {
    return(readAlleleCountFiles(paste0(Workdir,'/alleleCounts/',Worksheet$Patient_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'_unfiltered_chr'),'.txt',chrom_names,0))
  }
  names(NORMAL_COUNTS)=Worksheet$Normal_ID
  stopifnot(all(sapply(2:length(NORMAL_COUNTS), function(x) identical(rownames(NORMAL_COUNTS[[1]]),rownames(NORMAL_COUNTS[[x]])))))
  print('      Getting allelic information...')
  allele_data=readAllelesFiles(alleles.prefix,'.txt',chrom_names)
  stopifnot(identical(rownames(NORMAL_COUNTS[[1]]),rownames(allele_data)))
  # Get ref/alt/tot information for all cases based on alleleCounts and allele data
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) {
    x=x[,3:6]
    x$ref=x[cbind(1:nrow(allele_data),allele_data$a0)]
    x$alt=x[cbind(1:nrow(allele_data),allele_data$a1)]
    x$tot=x$ref+x$alt
    x$vaf=x$alt/x$tot
    x=x[,c('ref','alt','tot','vaf')]
    return(x)
  })
  # Get a DF with all total counts
  NORMAL_COUNTS_tot=do.call(cbind,lapply(NORMAL_COUNTS,function(x) x[,'tot',drop=F]))
  colnames(NORMAL_COUNTS_tot)=names(NORMAL_COUNTS)
  # Flag samples having too few SNPs (counts>minCounts)
  COVERED=apply(NORMAL_COUNTS_tot,2,function(x) length(which(x>=minCounts)))
  if (plotQC) {
    png(paste0(Workdir,'/plotQC/Low_number_of_SNPs.png'),height=15,width=15,units='cm',res=300,pointsize=6)
    par(mar=c(2.1,2.1,1.05,1.05))
    barplot(sort(COVERED),names.arg=NA,col='black',border=NA,space=0)
    abline(h=median(COVERED)/2,col='red')
    dev.off()
  }
  TO_REMOVE=names(which(COVERED<median(COVERED)/2))
  if (length(TO_REMOVE)>0) {
    print(paste0('   Remove samples with very low number of SNPs above threshold: ',paste(TO_REMOVE,collapse=', ')))
    NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[,setdiff(colnames(NORMAL_COUNTS_tot),TO_REMOVE)]
    NORMAL_COUNTS=NORMAL_COUNTS[colnames(NORMAL_COUNTS_tot)]
    Worksheet=Worksheet[-which(Worksheet$Normal_ID %in% TO_REMOVE),]
  }
  rm(TO_REMOVE,COVERED)
  # Flag samples having low coverages
  COVERED=apply(NORMAL_COUNTS_tot,2,function(x) sum(x[x>=minCounts]))
  if (plotQC) {
    png(paste0(Workdir,'/plotQC/Low_coverage.png'),height=15,width=15,units='cm',res=300,pointsize=6)
    par(mar=c(2.1,2.1,1.05,1.05))
    barplot(sort(COVERED),names.arg=NA,col='black',border=NA,space=0)
    abline(h=median(COVERED)/2,col='red')
    dev.off()
  }
  TO_REMOVE=names(which(COVERED<median(COVERED)/2))
  if (length(TO_REMOVE)>0) {
    print(paste0('   Remove samples with very low coverage: ',paste(TO_REMOVE,collapse=', ')))
    NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[,setdiff(colnames(NORMAL_COUNTS_tot),TO_REMOVE)]
    NORMAL_COUNTS=NORMAL_COUNTS[colnames(NORMAL_COUNTS_tot)]
    Worksheet=Worksheet[-which(Worksheet$Normal_ID %in% TO_REMOVE),]
  }
  rm(TO_REMOVE,COVERED)
  # Get SNPs covered in all samples
  IDs=rownames(NORMAL_COUNTS_tot[apply(NORMAL_COUNTS_tot,1,min)>0,])
  print(paste0('   Keep SNPs covered in all samples: ',length(IDs)))
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[IDs,])
  NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[IDs,]
  allele_data=allele_data[IDs,]
  rm(IDs)
  #######################################################################################
  # FILTER: keep SNPs with enough coverage (>=minCounts in at least 90% of cases) #
  #######################################################################################
  TO_KEEP=rowSums(NORMAL_COUNTS_tot>=minCounts)>=ncol(NORMAL_COUNTS_tot)*0.9
  print(paste0('   Keep SNPs with enough coverage in most cases: ',length(which(TO_KEEP)),' (-',round((nrow(NORMAL_COUNTS_tot)-length(which(TO_KEEP)))/nrow(NORMAL_COUNTS_tot)*100,2),'%)'))
  allele_data=allele_data[TO_KEEP,]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[TO_KEEP,])
  NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[TO_KEEP,]
  rm(TO_KEEP)
  # Get a DF with all VAFs
  NORMAL_COUNTS_vaf=do.call(cbind,lapply(NORMAL_COUNTS,function(x) x[,'vaf',drop=F]))
  colnames(NORMAL_COUNTS_vaf)=names(NORMAL_COUNTS)
  ############################################
  # FILTER: keep SNPs with finite BAF values #
  ############################################
  TO_KEEP=apply(NORMAL_COUNTS_vaf,1,function(x) all(is.finite(x)))
  print(paste0('   Keep SNPs with finite BAF values: ',length(which(TO_KEEP)),' (-',round((nrow(NORMAL_COUNTS_tot)-length(which(TO_KEEP)))/nrow(NORMAL_COUNTS_tot)*100,2),'%)'))
  allele_data=allele_data[TO_KEEP,]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[TO_KEEP,])
  NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[TO_KEEP,]
  NORMAL_COUNTS_vaf=NORMAL_COUNTS_vaf[TO_KEEP,]
  rm(TO_KEEP)
  #################################################################
  # FILTER: remove SNPs with noisy BAF values in multiple samples #
  #################################################################
  if (is.null(X_nonPAR) || length(which(Worksheet$Gender=='XX'))<10) {
    # No PAR/nonPAR information and/or too few females
    TO_KEEP=rowSums(NORMAL_COUNTS_tot>=minCounts&(NORMAL_COUNTS_vaf<0.9&NORMAL_COUNTS_vaf>0.68|NORMAL_COUNTS_vaf>0.1&NORMAL_COUNTS_vaf<0.32))<2.5*0.053876*rowSums(NORMAL_COUNTS_tot>=minCounts&NORMAL_COUNTS_vaf<0.9&NORMAL_COUNTS_vaf>0.1)
  } else {
    FEMALES=which(Worksheet$Gender=='XX')
    # Here, autosomes are 1:22 and PAR1/PAR2 regions in X. These regions should be 1+1
    INDEX_autosomes=which(allele_data$chromosome %in% 1:22 | (allele_data$chromosome=='X' & (allele_data$position<X_nonPAR[1] | allele_data$position>X_nonPAR[2])))
    # For autosomes, consider all samples since SNPs will be located in 1+1 regions
    TO_KEEP_autosomes=rowSums(NORMAL_COUNTS_tot[INDEX_autosomes,]>=minCounts&(NORMAL_COUNTS_vaf[INDEX_autosomes,]<0.9&NORMAL_COUNTS_vaf[INDEX_autosomes,]>0.68|NORMAL_COUNTS_vaf[INDEX_autosomes,]>0.1&NORMAL_COUNTS_vaf[INDEX_autosomes,]<0.32))<2.5*0.053876*rowSums(NORMAL_COUNTS_tot[INDEX_autosomes,]>=minCounts&NORMAL_COUNTS_vaf[INDEX_autosomes,]<0.9&NORMAL_COUNTS_vaf[INDEX_autosomes,]>0.1)
    # Here, X is only nonPAR (1+1 for females and 1+0 for males)
    INDEX_X=which(allele_data$chromosome=='X' & allele_data$position>=X_nonPAR[1] & allele_data$position<=X_nonPAR[2])
    if (length(INDEX_X)>0) {
      # For X, only consider females (1+1)
      TO_KEEP_X=rowSums(NORMAL_COUNTS_tot[INDEX_X,FEMALES]>=minCounts&(NORMAL_COUNTS_vaf[INDEX_X,FEMALES]<0.9&NORMAL_COUNTS_vaf[INDEX_X,FEMALES]>0.68|NORMAL_COUNTS_vaf[INDEX_X,FEMALES]>0.1&NORMAL_COUNTS_vaf[INDEX_X,FEMALES]<0.32))<2.5*0.053876*rowSums(NORMAL_COUNTS_tot[INDEX_X,FEMALES]>=minCounts&NORMAL_COUNTS_vaf[INDEX_X,FEMALES]<0.9&NORMAL_COUNTS_vaf[INDEX_X,FEMALES]>0.1)
      # Merge autosomes and X
      TO_KEEP=c(TO_KEEP_autosomes,TO_KEEP_X)[order(c(INDEX_autosomes,INDEX_X))]
      rm(TO_KEEP_X)
    } else {
      TO_KEEP=TO_KEEP_autosomes
    }
    rm(FEMALES,INDEX_autosomes,INDEX_X,TO_KEEP_autosomes)
  }
  
  print(paste0('   Remove SNPs with noisy BAF: ',length(which(TO_KEEP)),' (-',round((nrow(NORMAL_COUNTS_tot)-length(which(TO_KEEP)))/nrow(NORMAL_COUNTS_tot)*100,2),'%)'))
  allele_data=allele_data[TO_KEEP,]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[TO_KEEP,])
  NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[TO_KEEP,]
  NORMAL_COUNTS_vaf=NORMAL_COUNTS_vaf[TO_KEEP,]
  rm(TO_KEEP)
  ##################################################################
  # FILTER: remove SNPs close to each other with similar genotypes #
  ##################################################################
  SNP_dist=diff(allele_data[,2])
  TO_REMOVE = sapply(which(SNP_dist<=75),function(poske) {
    bafnp1 = NORMAL_COUNTS_vaf[poske,]
    bafnp2 = NORMAL_COUNTS_vaf[poske+1,]
    hetp1 = ifelse(bafnp1>0.1&bafnp1<0.9,1,0)
    hetp2 = ifelse(bafnp2>0.1&bafnp2<0.9,1,0)
    # if less than 10% of the calls are different
    if(sum(abs(hetp1-hetp2))<0.05*(sum(hetp1)+sum(hetp2))) {
      if(SNP_dist[max(poske-1,1)]<SNP_dist[min(poske+1,length(SNP_dist))]) {
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
    print(paste0('   Keep SNPs separated (distance and genotype): ',nrow(NORMAL_COUNTS_tot)-length(TO_REMOVE),' (-',round(length(TO_REMOVE)/nrow(NORMAL_COUNTS_tot)*100,2),'%)'))
    allele_data=allele_data[-TO_REMOVE,]
    NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[-TO_REMOVE,])
    NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[-TO_REMOVE,]
    NORMAL_COUNTS_vaf=NORMAL_COUNTS_vaf[-TO_REMOVE,]
  }
  rm(TO_REMOVE,SNP_dist)
  # Compute logR
  logR=computeLogR(NORMAL_COUNTS_tot)
  # Get standard deviation
  sdkes = apply(logR,1,sd)
  #######################################
  # FILTER: remove SNPs with noisy logR #
  #######################################
  TO_KEEP=sdkes<=2*median(sdkes)
  print(paste0('   Remove SNPs with noisy logR: ',length(which(TO_KEEP)),' (-',round((nrow(NORMAL_COUNTS_tot)-length(which(TO_KEEP)))/nrow(NORMAL_COUNTS_tot)*100,2),'%)'))
  allele_data=allele_data[TO_KEEP,]
  NORMAL_COUNTS=lapply(NORMAL_COUNTS,function(x) x[TO_KEEP,])
  NORMAL_COUNTS_tot=NORMAL_COUNTS_tot[TO_KEEP,]
  NORMAL_COUNTS_vaf=NORMAL_COUNTS_vaf[TO_KEEP,]
  rm(TO_KEEP,sdkes,logR)
  # Re-compute cleaned logR
  logR=computeLogR(NORMAL_COUNTS_tot)
  # Flag noisy sample(s)
  sdsam = apply(logR,2,sd)
  if (plotQC) {
    png(paste0(Workdir,'/plotQC/Noisy_logR.png'),height=15,width=15,units='cm',res=300,pointsize=6)
    par(mar=c(2.1,2.1,1.05,1.05))
    barplot(sort(sdsam),names.arg=NA,col='black',border=NA,space=0)
    abline(h=1.5*median(sdsam),col='red')
    dev.off()
  }
  IDs=names(which(sdsam>1.5*median(sdsam)))
  if (length(IDs)>0) {
    warning(paste0('Noisy sample',if (length(IDs)>1) 's',': ',paste0(IDs,collapse=', ')))
  }
  rm(IDs,sdsam)
  # Write & plot logR+BAF values
  if (!dir.exists(paste0(Workdir,'/Normal_data'))) dir.create(paste0(Workdir,'/Normal_data'))
  write.table(cbind(allele_data[,1:2],NORMAL_COUNTS_vaf),paste0(Workdir,"/Normal_data/Normal_BAF.txt"),sep="\t",col.names=NA,row.names=T,quote=F)
  write.table(cbind(allele_data[,1:2],logR),paste0(Workdir,"/Normal_data/Normal_LogR.txt"),sep="\t",col.names=NA,row.names=T,quote=F)
  plotLogRandBAF(paste0(Workdir,'/Normal_data'),allele_data[,1:2],logR,NORMAL_COUNTS_vaf)
  # Write output
  if (!dir.exists(paste0(Workdir,'/alleleData/Cleaned'))) dir.create(paste0(Workdir,'/alleleData/Cleaned'))
  for (CHR in chrom_names) {
    tmp=allele_data[which(allele_data[,1]==CHR),]
    write.table(tmp[,1:2],paste0(Workdir,'/alleleData/Cleaned/loci_chr',CHR,'.txt'),row.names=F,col.names=F,sep='\t',quote=F)
    write.table(tmp[,2:4],paste0(Workdir,'/alleleData/Cleaned/alleleData_chr',CHR,'.txt'),row.names=F,col.names=T,sep='\t',quote=F)
  }; rm(CHR,tmp)
}

#' Method to extract a curated list of SNPs covered by a targeted sequencing experiment. 
#'
#' From a complete set of loci (alleles.prefix), this method will keep SNPs falling into the targeted design (based on BED_file) and check allele counts in normal samples (listed in Worksheet). The cleaned list of loci/allele files will be located under Workdir/alleleData/Cleaned/.
#' 
#' @param Worksheet A table with the following columns: Patient_ID, Normal_ID, Normal_file and Gender. Must contain one single normal per patient. Normal_file can either be BAMs/CRAMs or paths to pre-computed (zipped) alleleCounts (e.g. "sample_alleleCounts_chr"). Gender must either be XX (females) or XY (males).
#' @param Workdir The folder where output should go (will be created if it doesn't exist).
#' @param alleles.prefix Prefix path to the allele data (e.g. "G1000_alleles_chr").
#' @param BED_file A BED file for only looking at SNPs within specific intervals. Must fit with the design used for targeted sequencing.
#' @param allelecounter_exe Path to the allele counter executable.
#' @param genomeVersion Genome version, either 'hg19' or 'hg38'.
#' @param nthreads The number of parallel processes to speed up the process (optional, default=1).
#' @param minCounts Minimum depth required in the normal for a SNP to be considered (optional, default=10).
#' @param chrom_names A vector containing the names of chromosomes to be considered (optional, default=c(1:22,'X')).
#' @param min_base_qual Minimum base quality required for a read to be counted (optional, default=20).
#' @param min_map_qual Minimum mapping quality required for a read to be counted (optional, default=35).
#' @param ref.fasta FASTA file used for generating CRAMs (optional, default=NA).
#' @param plotQC A boolean to generate QC reports as PNGs (optional, default=T).
#' @export
ascat.prepareTargetedSeq=function(Worksheet, Workdir, alleles.prefix, BED_file, allelecounter_exe, genomeVersion, nthreads=1,
                                  minCounts=10, chrom_names=c(1:22,'X'), min_base_qual=20, min_map_qual=35, ref.fasta=NA, plotQC=T) {
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  requireNamespace("foreach")
  requireNamespace("doParallel")
  
  ########
  # Init #
  ########
  stopifnot(genomeVersion %in% c('hg19','hg38'))
  if (genomeVersion=='hg19') {
    X_nonPAR=c(2699521,154931043)
  } else if (genomeVersion=='hg38') {
    X_nonPAR=c(2781480,155701382)
  }
  Worksheet=read.table(Worksheet,sep='\t',header=T,stringsAsFactors=F)
  stopifnot(nrow(Worksheet)>1)
  stopifnot(all(c('Patient_ID','Normal_ID','Normal_file','Gender') %in% colnames(Worksheet)))
  stopifnot(length(which(duplicated(Worksheet$Patient_ID)))==0)
  stopifnot(length(which(duplicated(Worksheet$Sample_ID)))==0)
  stopifnot(all(Worksheet$Gender %in% c('XX','XY')))
  registerDoParallel(cores=nthreads)
  if (!dir.exists(Workdir)) dir.create(Workdir)
  if (!dir.exists(paste0(Workdir,'/alleleData/Raw'))) dir.create(paste0(Workdir,'/alleleData/Raw'),recursive=T)
  if (!dir.exists(paste0(Workdir,'/alleleCounts'))) dir.create(paste0(Workdir,'/alleleCounts'))
  if (!dir.exists(paste0(Workdir,'/plotQC')) && plotQC) dir.create(paste0(Workdir,'/plotQC'))
  Process_HTS_file=F
  Process_AC_file=F
  
  # Two scenarios here: either Normal_file are BAMs/CRAMs or paths to pre-computed (zipped) alleleCounts
  if (all(sapply(strsplit(Worksheet$Normal_file,'\\.'),function(x) x[length(x)]) %in% c("bam",'BAM','cram','CRAM'))) {
    stopifnot(all(file.exists(Worksheet$Normal_file)))
    Process_HTS_file=T
  } else if (all(sapply(Worksheet$Normal_file,function(x) all(file.exists(paste0(x,chrom_names,'.txt')))))) {
    Process_AC_file=T
    SUFFIX='.txt'
  } else if (all(sapply(Worksheet$Normal_file,function(x) all(file.exists(paste0(x,chrom_names,'.txt.gz')))))) {
    Process_AC_file=T
    SUFFIX='.txt.gz'
  } else {
    stop('Worksheet$Normal_file seems incorrect; Please check that input files exist and have the right extension (bam/BAM/cram/CRAM from HTS data or txt/txt.gz for pre-computed allele counts.)')
  }
  
  # subset SNPs based on BED
  print('Subsetting SNPs based on BED')
  allele_data=readAllelesFiles(alleles.prefix,'.txt',chrom_names)
  BED=read.table(BED_file,sep='\t',header=F,stringsAsFactors=F)[,1:3]
  colnames(BED)=c('chr','start','end')
  BED$chr=gsub('^chr','',BED$chr)
  BED$start=BED$start+1 # Start is 0-based in BED files
  BED=BED[BED$chr %in% chrom_names,]
  if (nrow(BED)==0) stop('Major issue with BED file, please double-check its content')
  overlaps=GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames=BED$chr,ranges=IRanges::IRanges(start=BED$start,end=BED$end)),
                                       GenomicRanges::GRanges(seqnames=allele_data$chromosome,ranges=IRanges::IRanges(start=allele_data$position,end=allele_data$position)))
  allele_data=allele_data[unique(overlaps@to),]
  for (chr in chrom_names) {
    tmp=allele_data[allele_data$chromosome==chr,]
    tmp=tmp[order(tmp$position),]
    write.table(tmp[,-1],file=paste0(Workdir,'/alleleData/Raw/alleleData_chr',chr,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
    write.table(tmp[,1:2],file=paste0(Workdir,'/alleleData/Raw/loci_chr',chr,'.txt'),sep='\t',row.names=F,col.names=F,quote=F)
    rm(tmp)
  }; rm(chr,allele_data,BED,overlaps)
  
  if (Process_AC_file) {
    allele_data=readAllelesFiles(paste0(Workdir,'/alleleData/Raw/alleleData_chr'),'.txt',chrom_names)
    allele_data=lapply(chrom_names,function(x) rownames(allele_data[allele_data$chromosome==x,]))
    names(allele_data)=chrom_names
  }
  
  for (INDEX in 1:nrow(Worksheet)) {
    print(paste0('   Processing normal sample: ',Worksheet$Normal_ID[INDEX],' (',Worksheet$Patient_ID[INDEX],'; ',INDEX,'/',nrow(Worksheet),')'))
    if (!dir.exists(paste0(Workdir,'/alleleCounts/',Worksheet$Patient_ID[INDEX],'/',Worksheet$Normal_ID[INDEX]))) dir.create(paste0(Workdir,'/alleleCounts/',Worksheet$Patient_ID[INDEX],'/',Worksheet$Normal_ID[INDEX]),recursive=T)
    if (Process_HTS_file) {
      # Run alleleCounter on all samples (normals only)
      foreach(CHR=chrom_names) %dopar% {
        ascat.getAlleleCounts(seq.file=Worksheet$Normal_file[INDEX],
                              output.file=paste0(Workdir,'/alleleCounts/',Worksheet$Patient_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'_unfiltered_chr',CHR,'.txt'),
                              loci.file=paste0(Workdir,'/alleleData/Raw/loci_chr',CHR,'.txt'),
                              min.base.qual=min_base_qual,
                              min.map.qual=min_map_qual,
                              allelecounter.exe=allelecounter_exe,
                              ref.fasta=ref.fasta)
      }
    } else {
      foreach(CHR=chrom_names) %dopar% {
        LOCI=fread(paste0(Worksheet$Normal_file[INDEX],CHR,SUFFIX),sep='\t',data.table=F,showProgress=F,quote='')
        rownames(LOCI)=paste0(LOCI[,1],'_',LOCI[,2])
        stopifnot(all(allele_data[[as.character(CHR)]] %in% rownames(LOCI)))
        LOCI=LOCI[allele_data[[as.character(CHR)]],]
        write.table(LOCI,file=paste0(Workdir,'/alleleCounts/',Worksheet$Patient_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'/',Worksheet$Normal_ID[INDEX],'_unfiltered_chr',CHR,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
      }
    }
  }; rm(INDEX)
  
  ###########################
  # 3) Clean normal samples #
  ###########################
  print('Clean normal samples')
  getLociFromNormals(Worksheet=Worksheet,
                     Workdir=Workdir,
                     alleles.prefix=paste0(Workdir,'/alleleData/Raw/alleleData_chr'),
                     chrom_names=chrom_names,
                     minCounts=minCounts,
                     X_nonPAR=X_nonPAR,
                     plotQC=plotQC)
}