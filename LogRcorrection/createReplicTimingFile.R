suppressPackageStartupMessages(library(rtracklayer))

args = commandArgs(TRUE)
locifile = toString(args[1])

## load locifile for array and turn into GRanges object
loci <- read.delim(file = locifile, header = T, as.is = T)
loci_gr <- GRanges(seqnames = loci[,2], ranges = IRanges(start = loci[,3], width = 1))
names(loci_gr) <- loci[,1]
rm(loci)

## download and read wavelet smoothed repli-seq timing data from UCSC
UWrepliseqbigwigs <- c("wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig","wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig","wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig","wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig","wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig","wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig","wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig","wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig")
UCSCdownloadpath <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/"
print("Downloading UW Repli-seq data from UCSC")
invisible(lapply(X = UWrepliseqbigwigs, FUN = function(x) download.file(url = paste0(UCSCdownloadpath, x), destfile = x)))
wavelets <- lapply(UWrepliseqbigwigs, rtracklayer::import.bw)

## match each SNP locus with the nearest replication timing estimate from each cell line
seqlevelsStyle(wavelets[[1]]) <- "Ensembl"
nearestidxs <- nearest(x = loci_gr, subject = wavelets[[1]], select = "arbitrary")
mcols(loci_gr) <- do.call(cbind, lapply(X = wavelets, FUN = function(x) mcols(x)[nearestidxs, "score"]))
colnames(mcols(loci_gr)) <- sub(pattern = "^wgEncodeUwRepliSeq", replacement = "", x = sub(pattern = "WaveSignalRep1.bigWig$", replacement = "", x = UWrepliseqbigwigs))

locidf <- as.data.frame(loci_gr)[, -c(3:5)]
colnames(locidf) <- c("Chr", "Position", colnames(locidf)[-c(1,2)])
locidf[,3:ncol(locidf)]=round(locidf[,3:ncol(locidf)],6)
write.table(file = "ReplicationTiming_SNPloci.txt", x = locidf, sep = "\t", row.names = T, quote = F, col.names = NA)

## clean up
unlink(UWrepliseqbigwigs)
print("Done!")
