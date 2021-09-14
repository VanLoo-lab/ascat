process=function(CHRSEQ,SNP_POS,WINDOW,THRESH,NCORES) {
  CHRSIZE=length(CHRSEQ)
  # Quick check: make sure that A/C/G/T/N represent >99.99% of reference sequence
  stopifnot(letterFrequency(CHRSEQ,c('ACGTN'))/CHRSIZE>0.9999)
  if (WINDOW%%2==0) WINDOW=WINDOW+1
  # Add Ns at chromosome extremities so nearby SNPs can be processed.
  CHRSEQ=xscat(paste0(rep('N',floor(WINDOW/2)),collapse=''),CHRSEQ,paste0(rep('N',floor(WINDOW/2)),collapse=''))
  # Get G+C and A+T counts for all positions
  LFISV=as.data.frame(letterFrequencyInSlidingView(CHRSEQ,WINDOW,c('GC','AT'),OR='_'))
  # Make sure that we have one value for each position
  stopifnot(CHRSIZE==nrow(LFISV))
  # Keep SNP information
  LFISV=LFISV[SNP_POS$Position,]
  # Effective window size (A/C/G/T only)
  LFISV$NREALS=LFISV$G_C+LFISV$A_T
  # Compute GC%
  LFISV$GC=LFISV$G_C/LFISV$NREALS
  # If effective size is smaller than threshold, set NA
  LFISV$GC[LFISV$NREALS<THRESH]=NA
  return(round(LFISV$GC,6))
}

suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(Biostrings))
options(scipen=999,digits=6)
THRESH=20 # Defines the minimum number of A/C/G/T to compute GC%, otherwise set NA
WINDOWS=c(25,50,100,200,500,1e3,2e3,5e3,1e4,
          2e4,5e4,1e5,2e5,5e5,1e6)
names(WINDOWS)=c('25bp','50bp','100bp','200bp','500bp','1kb','2kb','5kb','10kb',
                 '20kb','50kb','100kb','200kb','500kb','1Mb')

args=commandArgs(trailingOnly=TRUE)
LOCI_file=args[1]
stopifnot(file.exists(LOCI_file))
LOCI=read.table(LOCI_file,sep='\t',header=T,stringsAsFactors=F,row.names=1)
NCORES=as.numeric(args[2])
stopifnot(is.finite(NCORES) && NCORES %in% 1:24)
registerDoParallel(cores=NCORES)
refFASTA=args[3]
stopifnot(file.exists(refFASTA))
print('Loading reference genome')
CHRSEQ=readDNAStringSet(refFASTA,format='fasta')
names(CHRSEQ)=sapply(strsplit(names(CHRSEQ),' '),function(x) x[[1]])
CHRSEQ=CHRSEQ[which(names(CHRSEQ) %in% as.character(LOCI$Chr))]
stopifnot(all(unique(LOCI$Chr) %in% names(CHRSEQ)))

OUT=foreach(CHR=unique(as.character(LOCI$Chr)),.combine=rbind) %dopar% {
  print(paste0('Processing chromosome: ',CHR))
  stopifnot(CHR %in% names(CHRSEQ))
  SNP_POS=LOCI[which(LOCI$Chr==CHR),]
  for (i in 1:length(WINDOWS)) {
    SNP_POS[,names(WINDOWS)[i]]=process(CHRSEQ[[CHR]],SNP_POS,WINDOWS[i],THRESH,NCORES)
  }; rm(i)
  return(SNP_POS)
}
write.table(OUT,file='GCcontent_SNPloci.txt',col.names=NA,row.names=T,quote=F,sep='\t')
print("Done!")