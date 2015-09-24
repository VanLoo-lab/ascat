SNPfile=args[4]
bedgraphs=args[5]

SNPpos = read.table(SNPfile, sep="\t",header=T,row.names=1)
platform = gsub(".txt","",SNPfile)

for (chr in c(1:22,"X")) {
  chrpos = SNPpos[SNPpos[,1]==chr,2]
  pattern = paste(chrpos,"\t",sep="")
  write.table(pattern,paste(platform,".chr",chr,".txt",sep=""),sep="\t", col.names=F, row.names=F,quote=F)
}


window1 = c("12","25","50",
            "100","250","500",
            "1000","2500","5000",
            "10000","25000","50000",
            "100000","250000","500000",
            "1000000","2500000","5000000")
window2 = c("25","50","100",
            "200","500","1000",
            "2000","5000","10000",
            "20000","50000","100000",
            "200000","500000","1000000",
            "2000000","5000000","10000000")
window3 = c("25bp","50bp","100bp",
            "200bp","500bp","1kb",
            "2kb","5kb","10kb",
            "20kb","50kb","100kb",
            "200kb","500kb","1Mb",
            "2Mb","5Mb","10Mb")
# memlimit = c("","","",
#              "","","",
#              "","","",
#              "","","",
#              
#              "","-R\"select[mem>500] rusage[mem=500]\" -M500","-R\"select[mem>500] rusage[mem=500]\" -M500",
#              "-R\"select[mem>500] rusage[mem=500]\" -M500","-R\"select[mem>500] rusage[mem=500]\" -M500",
#                 "-R\"select[mem>1000] rusage[mem=1000]\" -M1000")  

torunall = NULL
for (chr in c(1:22,"X")) {
#   torun = paste("bsub -q normal ",memlimit,"-J \"GC.chr",chr,".",window3,"\" -o GC.chr",chr,".",window3,
#                 ".out.txt \"perl make_windowed_feature_means.pl -w ",window1," -chr chr",chr,
#                 " -f /srv/data/vanloo/download_documentation/processed_bedgraphs/chr",
#                 chr,"/hg19.gc5Base.chr",chr,".bg.gz | fgrep -f ",platform,".chr",chr,".txt > GCcontent.chr",chr,
#                 ".",window3,".txt\"",sep="")
  
  torun = paste("perl make_windowed_feature_means.pl -w ",window1," -chr chr",chr,
                " -f ",bedgraphs,"/hg19.gc5Base.chr",chr,".bg.gz | fgrep -f ",platform,".chr",chr,".txt > GCcontent.chr",chr,
                ".",window3,".txt",sep="")
  torunall = c(torunall,torun)
}
write.table(torunall,"torun.sh",sep="\t",col.names=F,row.names=F,quote=F)

system("./torun.sh")

# there is a bug in the grep statement in that more lines can get returned than needed.. 
# This code should be robust to that though..
# one way to solve this is to add the chromosome!

# build matrix
GCcontent = matrix(nrow = dim(SNPpos)[1], ncol = length(window3))
rownames(GCcontent) = rownames(SNPpos)
colnames(GCcontent) = window3
for (chr in c(1:22,"X")) {
  GCcontentchr = GCcontent[SNPpos[,1]==chr,]
  posmapping = rownames(SNPpos)[SNPpos[,1]==chr]
  names(posmapping) = SNPpos[SNPpos[,1]==chr,2]
  for (i in 1:length(window3)) {
    GC = read.table(paste("GCcontent.chr",chr,".",window3[i],".txt",sep=""),header=F,row.names=1,sep="\t")
    GCcontentchr[posmapping,i] = GC[names(posmapping),1] 
  }
  GCcontent[rownames(GCcontentchr),] = GCcontentchr
}

write.table(cbind(SNPpos,GCcontent), paste(platform,".txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
