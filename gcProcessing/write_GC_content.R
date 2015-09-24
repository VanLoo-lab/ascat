#setwd("/lustre/scratch104/sanger/tjm/GC_correct_files")


args = commandArgs(T)
indices = as.integer(args[1])
g1000_prefix = toString(args[2])
gc_perc_prefix = toString(args[3])
outfile_prefix = toString(args[4])
if (indices==23) {indices = "X"}

chr_num <- c(1:23)
window3 = c("25bp","50bp","100bp","200bp","500bp","1kb","2kb","5kb","10kb","20kb","50kb","100kb","200kb","500kb","1Mb","2Mb","5Mb","10Mb")

# build matrix
for (chr in indices) {
  #SNPpos = read.table(paste("/lustre/scratch112/sanger/cgppipe/drj_tmp_caveman_test_08_08_2014/tmp/1000genomesloci/1000genomesloci2012_chr",chr_num[c(1:22,"X")%in%chr],".txt", sep=""), sep="\t")
	print(paste(g1000_prefix,chr_num[c(1:22,"X")%in%chr],".txt", sep=""))
  SNPpos = read.table(paste(g1000_prefix,chr_num[c(1:22,"X")%in%chr],".txt", sep=""), sep=" ")
  colnames(SNPpos) <- c("chr", "Position")
  GCcontent = matrix(nrow = dim(SNPpos)[1], ncol = length(window3))  
  rownames(GCcontent) = SNPpos[,2]
  colnames(GCcontent) = window3
  for (i in 1:length(window3)) {
    print(paste("starting window",i))
    #GC = read.table(paste("/lustre/scratch104/sanger/tjm/GC_correct_files/GCcontent.chr",chr,".",window3[i],".txt",sep=""),header=F,sep="\t")
  print(paste(gc_perc_prefix,chr,".",window3[i],".txt",sep=""))
    GC = read.table(paste(gc_perc_prefix,chr,".",window3[i],".txt",sep=""),header=F,sep="\t")
    print(paste("have read GC window",i))
    GCcontent[which(SNPpos[,2]%in%GC[,3]),i] = GC[which(GC[,3]%in%SNPpos[,2]),4] 
  }
  print("about to write output")
  #write.table(cbind(SNPpos,GCcontent),paste("/lustre/scratch104/sanger/tjm/GC_correct_files/1000_genomes_GC_corr_chr_",chr,".txt",sep=""),sep="\t",quote=F)
  write.table(cbind(SNPpos,GCcontent),paste(outfile_prefix,chr,".txt",sep=""),sep="\t",quote=F)
}

