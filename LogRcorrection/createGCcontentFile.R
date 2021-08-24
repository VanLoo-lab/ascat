process=function(CHRSEQ,SNP_POS,WINDOW,THRESH,NCORES) {
  subprocess=function(CHRSEQ,START,STOP,SIZE,THRESH) {
    VALS=letterFrequency(DNAStringSet(CHRSEQ,start=START,end=STOP),c('G','C','N'))
    NREALS=SIZE-VALS[,'N']
    GC=(VALS[,'G']+VALS[,'C'])/NREALS
    GC[NREALS<THRESH]=NA
    return(GC)
  }
  FLANK=ceiling((WINDOW-1)/2)
  SIZE=FLANK*2+1
  SNP_POS$START=SNP_POS$Position-FLANK
  SNP_POS$STOP=SNP_POS$Position+FLANK
  INDEX=which(SNP_POS$START<1); if (length(INDEX)>0) {SNP_POS$START[INDEX]=1}; rm(INDEX)
  INDEX=which(SNP_POS$STOP>CHRSEQ@length); if (length(INDEX)>0) {SNP_POS$STOP[INDEX]=CHRSEQ@length}; rm(INDEX)
  SEQS=round(seq(1,nrow(SNP_POS),length.out=NCORES+1))
  SEQS=lapply(1:NCORES,function(x) c(SEQS[x],ifelse(x==NCORES,SEQS[x+1],SEQS[x+1]-1)))
  GC=foreach(x=1:NCORES) %dopar%  subprocess(CHRSEQ,SNP_POS$START[(SEQS[[x]][1]):(SEQS[[x]][2])],SNP_POS$STOP[(SEQS[[x]][1]):(SEQS[[x]][2])],SIZE,THRESH)
  rm(SNP_POS)
  gc()
  return(round(as.numeric(unlist(GC)),6))
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
CHRSEQ=readDNAStringSet(refFASTA,format='fasta')
names(CHRSEQ)=sapply(strsplit(names(CHRSEQ),' '),function(x) x[[1]])
CHRSEQ=CHRSEQ[which(names(CHRSEQ) %in% as.character(LOCI$Chr))]
stopifnot(all(unique(LOCI$Chr) %in% names(CHRSEQ)))

OUT=foreach(CHR=unique(as.character(LOCI$Chr)),.combine=rbind) %do% {
  print(paste0('Processing chromosome: ',CHR))
  stopifnot(CHR %in% names(CHRSEQ))
  SNP_POS=LOCI[which(LOCI$Chr==CHR),]
  for (i in 1:length(WINDOWS)) {
    print(paste0('   Window: ',names(WINDOWS)[i]))
    SNP_POS[,names(WINDOWS)[i]]=process(CHRSEQ[[CHR]],SNP_POS,WINDOWS[i],THRESH,NCORES)
  }; rm(i)
  return(SNP_POS)
}
write.table(OUT,file='GCcontent_SNPloci.txt',col.names=NA,row.names=T,quote=F,sep='\t')
print("Done!")