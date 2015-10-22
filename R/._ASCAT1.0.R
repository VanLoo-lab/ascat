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
  psi_pos = seq(1,5.4,0.05) 
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
      d[i,j] = sum(abs(nA - pmax(round(nA),0))^2 + abs(nB - pmax(round(nB),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T) 
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
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, distancepng = NA, copynumberprofilespng = NA, gamma = 0.55) {
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

  axis(1, at = seq(0, 4/4.4, by = 1/4.4), label = seq(1, 5, by = 1))
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

        goodnessOfFit = (1-m/TheoretMaxdist) * 100

        if (ploidy > 1.6 & ploidy < 4.8 & rho >= 0.2 & goodnessOfFit > 80 & (percentzero > 0.01 | perczeroAbb > 0.1)) {
          nropt = nropt + 1
          optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
          localmin[nropt] = m
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
        ploidy_opt1 = optima[[i]][4]
        goodnessOfFit_opt1 = optima[[i]][5]
        points((psi_opt1-1)/4.4,(rho_opt1-0.1)/0.95,col="green",pch="X", cex = 2)
      }
    }
  }

  if (!is.na(distancepng)) {
    dev.off()
  }


  if(nropt>0) {

    if (!is.na(copynumberprofilespng)) {
      png(filename = copynumberprofilespng, width = 2000, height = 1000, res = 200)
    } 
    else {      
      windows(10,5)
    }

    par(mar = c(0.5,5,5,0.5), mfrow=c(2,1), cex = 0.4, cex.main=3, cex.axis = 2.5)

    rho = rho_opt1
    psi = psi_opt1
    
    nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
    plot(c(1001,length(nAfull)-1000), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nA-0.1,col="red",pch = "|")
    points(nB+0.1,col="green",pch = "|")
# don't ask me why, but the "|" sign is not centered, so the lines need to be shifted..
    abline(v=-20,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos-20,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
      abline(v=vpos-20,lty=1,col="lightgrey")
    }
    
    rBacktransform = gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
    bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
    rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))
    maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
    plot(c(1001,length(nAfull)-1000), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(confidence,col="blue",pch = "|")
    abline(v=-20,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos-20,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
      abline(v=vpos-20,lty=1,col="lightgrey")
    }
    

    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }

    rho = rho_opt1
    psi = psi_opt1

    nA = nAfull
    nB = nBfull

    probesWithoutNA = intersect(names(as.matrix(na.omit(r))[,]),names(as.matrix(na.omit(b))[,]))
    bBeforePCF = baf[probesWithoutNA]
    names(bBeforePCF) = probesWithoutNA

    n1hetero = pmax(ifelse(bBeforePCF<0.5,round(nA),round(nB)),0)
    n2hetero = pmax(ifelse(bBeforePCF>=0.5,round(nA),round(nB)),0)

    n1all = vector(mode = "integer", length = length(lrr))
    n1all[] = NA
    names(n1all) = names(lrr)
    n1all[probesWithoutNA] = n1hetero

    n2all = vector(mode = "integer", length = length(lrr))
    n2all[] = NA
    names(n2all) = names(lrr)
    n2all[probesWithoutNA] = n2hetero

    homoProbesWithoutNA = setdiff(names(as.matrix(na.omit(lrrsegmented))[,]),names(bafsegmented))
    homoBBeforePCF = baf[homoProbesWithoutNA]

    rhomo = lrrsegmented[homoProbesWithoutNA]

    nA = (2*rho-2+2^(rhomo/gamma)*((1-rho)*2+rho*psi))/rho
    nB = rep(0,length(nA))

    n1homo = pmax(ifelse(homoBBeforePCF<0.5,round(nA),round(nB)),0)
    n2homo = pmax(ifelse(homoBBeforePCF>=0.5,round(nA),round(nB)),0)

    n1all[homoProbesWithoutNA] = n1homo
    n2all[homoProbesWithoutNA] = n2homo

    return(list(rho = rho_opt1, psi = psi_opt1, nA = n1all, nB = n2all))
  }

  else {
    return(list(rho = NA, psi = NA, nA = NA, nB = NA))
  }
  
}


