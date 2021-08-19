#' @title ascat.GCcorrect
#' @description Corrects logR of the tumour sample(s) with genomic GC content (replication timing is optional)
#' @param ASCATobj an ASCAT object
#' @param GCcontentfile File containing the GC content around every SNP for increasing window sizes
#' @param replictimingfile File containing replication timing at every SNP for various cell lines (optional)
#' @details Note that probes not present in the GC content file will be lost from the results
#' @return ASCAT object with corrected tumour logR
#'
#' @export
ascat.GCcorrect = function(ASCATobj, GCcontentfile = NULL, replictimingfile = NULL) {
  if (is.null(GCcontentfile)) {
    stop("No GC content file given!")
  } else {
    GC_newlist=readCorrectionFile(GCcontentfile)
    stopifnot(!is.null(GC_newlist))
    ovl=intersect(rownames(ASCATobj$Tumor_LogR),rownames(GC_newlist))
    stopifnot(length(ovl)>nrow(ASCATobj$Tumor_LogR)/10)
    GC_newlist=GC_newlist[ovl,]
    if (!is.null(replictimingfile)) {
      replic_newlist=readCorrectionFile(replictimingfile)
      stopifnot(!is.null(replic_newlist))
      stopifnot(all(ovl %in% rownames(replic_newlist)))
      replic_newlist=replic_newlist[ovl,]
    } else {
      print.noquote("Warning: no replication timing file given, proceeding with GC correction only!")
    }
    
    SNPpos = ASCATobj$SNPpos[ovl,]
    Tumor_LogR = ASCATobj$Tumor_LogR[ovl,,drop=F]
    Tumor_BAF = ASCATobj$Tumor_BAF[ovl,,drop=F]
    
    chrs = intersect(ASCATobj$chrs,unique(SNPpos[,1]))
    
    Germline_LogR = NULL
    Germline_BAF = NULL
    if(!is.null(ASCATobj$Germline_LogR)) {
      Germline_LogR = ASCATobj$Germline_LogR[ovl,,drop=F]
      Germline_BAF = ASCATobj$Germline_BAF[ovl,,drop=F]
    }
    
    last = 0;
    ch = list();
    for (i in 1:length(ASCATobj$chrs)) {
      chrke = SNPpos[SNPpos[,1]==ASCATobj$chrs[i],]
      chrpos = chrke[,2]
      names(chrpos) = rownames(chrke)
      chrpos = sort(chrpos)
      ch[[i]] = (last+1):(last+length(chrpos))
      last = last+length(chrpos)
    }
    
    GC_correction_before=c()
    GC_correction_after=c()
    RT_correction_before=c()
    RT_correction_after=c()
    
    for (s in 1:length(ASCATobj$samples)) {
      print.noquote(paste("Sample ", ASCATobj$samples[s], " (",s,"/",length(ASCATobj$samples),")",sep=""))
      Tumordata = Tumor_LogR[,s]
      names(Tumordata) = rownames(Tumor_LogR)
      
      # Calculate correlation (explicit weighting now automatically done)
      corr = abs(cor(GC_newlist[, 3:ncol(GC_newlist)], Tumordata, use="complete.obs")[,1])
      
      index_1kb = grep(pattern = "X?1([06]00bp|kb)", x = names(corr))
      maxGCcol_insert = names(which.max(corr[1:index_1kb]))
      # if no replication timing data, expand large windows to 500kb
      if (!is.null(replictimingfile)) {
        index_max = grep(pattern = "X?10((00|24)00bp|0kb)", x = names(corr))
      } else {
        index_max = grep(pattern = "X?1M", x = names(corr)) - 1
      }
      # start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
      maxGCcol_amplic = names(which.max(corr[(index_1kb+2):index_max]))
      
      cat("GC correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
      cat("Short window size: ",maxGCcol_insert,"\n")
      cat("Long window size: ",maxGCcol_amplic,"\n")
      
      corrdata = data.frame(logr = Tumordata,
                            GC_insert = GC_newlist[,maxGCcol_insert],
                            GC_amplic = GC_newlist[,maxGCcol_amplic])
      
      if (!is.null(replictimingfile)) {
        corr_rep = abs(cor(replic_newlist[, 3:ncol(replic_newlist)], Tumordata, use="complete.obs")[,1])
        maxreplic = names(which.max(corr_rep))
        
        cat("Replication timing correlation: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n") 
        cat("Replication dataset: " ,maxreplic,"\n")
        
        # Multiple regression 
        corrdata$replic <- replic_newlist[, maxreplic]
        model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T) + splines::ns(x = replic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
        Tumor_LogR[,s] = residuals(model)
        
        RT_correction_before=c(RT_correction_before,paste0(maxreplic,"=",round(corr_rep[maxreplic],4)))
        corr_rep_after = abs(cor(replic_newlist[, which(colnames(replic_newlist)==maxreplic),drop=F], Tumor_LogR[,s], use="complete.obs")[,1])
        RT_correction_after=c(RT_correction_after,paste0(maxreplic,"=",round(corr_rep_after[maxreplic],4)))
      } else {
        model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
        Tumor_LogR[,s] = residuals(model)
        RT_correction_before=c(RT_correction_before,NA)
        RT_correction_after=c(RT_correction_after,NA)
      }
      
      chr = ASCAT:::split_genome(SNPpos)
      
      GC_correction_before=c(GC_correction_before,paste0(maxGCcol_insert,"=",round(corr[maxGCcol_insert],4),' / ',maxGCcol_amplic,'=',round(corr[maxGCcol_amplic],4)))
      corr_after=abs(cor(GC_newlist[, which(colnames(GC_newlist) %in% c(maxGCcol_insert,maxGCcol_amplic))], Tumor_LogR[,s], use="complete.obs")[,1])
      GC_correction_after=c(GC_correction_after,paste0(maxGCcol_insert,"=",round(corr_after[maxGCcol_insert],4),' / ',maxGCcol_amplic,'=',round(corr_after[maxGCcol_amplic],4)))
    }
    
    names(GC_correction_before)=colnames(Tumor_LogR)
    names(GC_correction_after)=colnames(Tumor_LogR)
    names(RT_correction_before)=colnames(Tumor_LogR)
    names(RT_correction_after)=colnames(Tumor_LogR)
    return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
                Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
                Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
                SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
                samples = colnames(Tumor_LogR), gender = ASCATobj$gender,
                sexchromosomes = ASCATobj$sexchromosomes,
                GC_correction_before=GC_correction_before, GC_correction_after=GC_correction_after,
                RT_correction_before=RT_correction_before, RT_correction_after=RT_correction_after))
  }
}

#' Function to read any type of correction file (should have very similar format: SNP ID, Chr, Position and then data)
#' @noRd
readCorrectionFile=function(correction_file) {
  if (file.exists(correction_file) && file.info(correction_file)$size>0) {
    header=as.character(readr::read_delim(file=correction_file, delim='\t', col_names=F, n_max=1, skip=0, col_types = readr::cols()))[-1]
    tmp=readr::read_tsv(file=correction_file,col_names=F,progress=F,col_types=paste0('cci',paste(rep('d',length(header)-2),collapse='')),skip=1)
    tmp=data.frame(tmp,stringsAsFactors=F)
    rownames(tmp)=tmp[,1]
    tmp[,1]=NULL
    colnames(tmp)=header
    colnames(tmp)[1:2]=c("Chr","Position")
    return(tmp)
  } else {
    return(NULL)
  }
}