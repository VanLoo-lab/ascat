# 
# Enhanced bivariate PCF filter for aCGH data (v. 08.02.2010)
# Whole chromosomes/chromosome arms wrapper function
#

fastAspcf <- function(logR, allB, kmin, gamma){

  N <- length(logR)
  w <- 1000 #w: windowsize
  d <- 100

  startw = -d
  stopw = w-d
  
  nseg = 0
  var2 = 0
  breakpts = 0
  larger = TRUE
  repeat{
    from <- max(c(1,startw))
    to  <-  min(c(stopw,N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]

    sd1 <- getMad(logRpart)
    sd2 <- getMad(allBflip)

    #Must check that sd1 and sd2 are defined and != 0:
    sd.valid <- c(!is.na(sd1),!is.na(sd2),sd1!=0,sd2!=0)
    if(all(sd.valid)){
      #run aspcfpart:
      #part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      breakptspart <- part.res$breakpts
      # the 'larger' is (occasionally) necessary in the last window of the segmentation!
      larger = breakptspart>breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      nseg = nseg+1
    }
    
    if(stopw < N+d){
      startw <- min(stopw-2*d + 1,N-2*d)
      stopw <- startw + w
    }else{
      break
    }
    
  }#end repeat
  breakpts <- unique(c(breakpts, N))
  if(nseg==0){nseg=1}  #just in case the sd-test never passes.
  sd2 <- sqrt(var2/nseg)
  
  # On each segment calculate mean of unflipped B allele data
  frst <- breakpts[1:length(breakpts)-1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)
  	
	yhat1 <- rep(NA,N)
  yhat2 <- rep(NA,N)

  for(i in 1:nseg){ 
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i]-frst[i]+1)
    yi2 <- allB[frst[i]:last[i]]
    # Center data around zero (by subtracting 0.5) and estimate mean
    if(length(yi2)== 0){
      mu <- 0
    }else{
      mu <- mean(abs(yi2-0.5))
    }
    
    # Make a (slightly arbitrary) decision concerning branches
    # This may be improved by a test of equal variances
    if(sqrt(sd2^2+mu^2) < 2*sd2){
      mu <- 0
    }
    yhat2[frst[i]:last[i]] <- rep(mu+0.5,last[i]-frst[i]+1)
  }

  return(list(yhat1=yhat1,yhat2=yhat2))

}#end fastAspcf



aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma){

	from <- max(c(1,a))
  usefrom <- max(c(1,a+d))
  useto <- min(c(N,b-d))
  
  N <- length(logRpart)
	y1 <- logRpart	
	y2 <- allBflip
	
	#Check that vectors are long enough to run algorithm:
	if(N < 2*kmin){
	 breakpts <- 0
   return(list(breakpts=breakpts))  
	}
 
	# Find initSum, initKvad, initAve for segment y[1..kmin]
	initSum1 <- sum(y1[1:kmin])
	initKvad1 <- sum(y1[1:kmin]^2)
	initAve1 <- initSum1/kmin     
	initSum2 <- sum(y2[1:kmin])
	initKvad2 <- sum(y2[1:kmin]^2)
	initAve2 <- initSum2/kmin

	# Define vector of best costs
	bestCost <- rep(0,N)
	cost1 <- (initKvad1 - initSum1*initAve1)/sd1^2
	cost2 <- (initKvad2 - initSum2*initAve2)/sd2^2
	bestCost[kmin] <- cost1 + cost2

	# Define vector of best splits
	bestSplit <- rep(0,N)
		
	# Define vector of best averages
	bestAver1 <- rep(0,N)
	bestAver2 <- rep(0,N)
	bestAver1[kmin] <- initAve1
	bestAver2[kmin] <- initAve2

	
	#Initialize
	Sum1 <- rep(0,N)
  Sum2 <- rep(0,N)
	Kvad1 <- rep(0,N)
  Kvad2 <- rep(0,N)
  Aver1 <- rep(0,N)
  Aver2 <- rep(0,N)
  Cost <- rep(0,N)
  
  # We have to treat the region y(1..2*kmin-1) separately, as it
	# cannot be split into two full segments
  kminP1 <- kmin+1
  for (k in (kminP1):(2*kmin-1)) {
		Sum1[kminP1:k] <- Sum1[kminP1:k]+y1[k]
		Aver1[kminP1:k] <- Sum1[kminP1:k]/((k-kmin):1)
		Kvad1[kminP1:k] <- Kvad1[kminP1:k]+y1[k]^2
		Sum2[kminP1:k] <- Sum2[kminP1:k]+y2[k]
		Aver2[kminP1:k] <- Sum2[kminP1:k]/((k-kmin):1)
		Kvad2[kminP1:k] <- Kvad2[kminP1:k]+y2[k]^2
		
		
		bestAver1[k] <- (initSum1+Sum1[kminP1])/k
		bestAver2[k] <- (initSum2+Sum2[kminP1])/k
    cost1 <- ((initKvad1+Kvad1[kminP1])-k*bestAver1[k]^2)/sd1^2
    cost2 <- ((initKvad2+Kvad2[kminP1])-k*bestAver2[k]^2)/sd2^2
  
    bestCost[k] <- cost1 + cost2
	  
 	}
  

  for (n in (2*kmin):N) {
		
   		nMkminP1=n-kmin+1
   		
   		Sum1[kminP1:n] <- Sum1[kminP1:n]+ y1[n]
   		Aver1[kminP1:n] <- Sum1[kminP1:n]/((n-kmin):1)
   		Kvad1[kminP1:n] <- Kvad1[kminP1:n]+ (y1[n])^2
   		                          
   		cost1 <- (Kvad1[kminP1:nMkminP1]-Sum1[kminP1:nMkminP1]*Aver1[kminP1:nMkminP1])/sd1^2
   		
      Sum2[kminP1:n] <- Sum2[kminP1:n]+ y2[n]
   		Aver2[kminP1:n] <- Sum2[kminP1:n]/((n-kmin):1)
   		Kvad2[kminP1:n] <- Kvad2[kminP1:n]+ (y2[n])^2
   		cost2 <- (Kvad2[kminP1:nMkminP1]-Sum2[kminP1:nMkminP1]*Aver2[kminP1:nMkminP1])/sd2^2
   		
      Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2
      
      Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
   		cost <- Cost[Pos] + gamma
   		
   		aver1 <- Aver1[Pos]
   		aver2 <- Aver2[Pos]
   		totAver1 <- (Sum1[kminP1]+initSum1)/n
   		totCost1 <- ((Kvad1[kminP1]+initKvad1) - n*totAver1*totAver1)/sd1^2
      totAver2 <- (Sum2[kminP1]+initSum2)/n
      totCost2 <- ((Kvad2[kminP1]+initKvad2) - n*totAver2*totAver2)/sd2^2
      totCost <- totCost1 + totCost2
   		
      if (totCost < cost) {
       		Pos <- 1
          cost <- totCost
          aver1 <- totAver1
          aver2 <- totAver2
   		}
   		bestCost[n] <- cost
   		bestAver1[n] <- aver1
   		bestAver2[n] <- aver2
   		bestSplit[n] <- Pos-1
 	
	  
   }#endfor
 	

	# Trace back
	n <- N
	breakpts <- n
	while(n > 0){
	  breakpts <- c(bestSplit[n], breakpts)
		n <- bestSplit[n]
	}#endwhile

  breakpts <- breakpts + from -1
  breakpts <- breakpts[breakpts>=usefrom & breakpts<=useto] 
  
  return(list(breakpts=breakpts))
  
}#end aspcfpart






div <- function(a, b, c){

	if(nargs() < 3){
	 	c <- 0
	 }#endif
 
	if(b > 0){
		v <-  a/b
	}else{
    v <- c
	}#endif
 
	 return(v)
}#endfunction


#function v = div(a, b, c)
# if nargin < 3
#     c = 0;
# end
# if b > 0
#     v = a/b;
# else
#     v = c;
# end
#end



#Get mad SD (based on KL code)
getMad <- function(x,k=25){
  
  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]
  
  #Calculate runMedian  
  runMedian <- medianFilter(x,k)
  
  dif <- x-runMedian
  SD <- mad(dif)
 
	return(SD)
}





exactPcf <- function(y, kmin, gamma) {
## Implementaion of exact PCF by Potts-filtering
	## x: input array of (log2) copy numbers
	## kmin: Mininal length of plateaus
	## gamma: penalty for each discontinuity
  	N <- length(y)
  	yhat <- rep(0,N);
  	if (N < 2*kmin) {
    		yhat <- rep(mean(y),N)
    		return(yhat)
  	}
  	initSum <- sum(y[1:kmin])
  	initKvad <- sum(y[1:kmin]^2)
  	initAve <- initSum/kmin;
  	bestCost <- rep(0,N)
  	bestCost[kmin] <- initKvad - initSum*initAve
  	bestSplit <- rep(0,N)
  	bestAver <- rep(0,N)
  	bestAver[kmin] <- initAve
  	Sum <- rep(0,N)
  	Kvad <- rep(0,N)
  	Aver <- rep(0,N)
  	Cost <- rep(0,N)
  	kminP1=kmin+1
  	for (k in (kminP1):(2*kmin-1)) {
    		Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    		Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
    		Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
    		bestAver[k] <- (initSum+Sum[kminP1])/k
    		bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
  	}
  	for (n in (2*kmin):N) {
   		yn <- y[n]
   		yn2 <- yn^2
   		Sum[kminP1:n] <- Sum[kminP1:n]+yn
   		Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
   		Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
   		nMkminP1=n-kmin+1
   		Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
   		Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
   		cost <- Cost[Pos]
   		aver <- Aver[Pos]
   		totAver <- (Sum[kminP1]+initSum)/n
   		totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
   		if (totCost < cost) {
       		Pos <- 1
          cost <- totCost
          aver <- totAver
   		}
   		bestCost[n] <- cost
   		bestAver[n] <- aver
   		bestSplit[n] <- Pos-1
 	}
 	n <- N
 	while (n > 0) {
   		yhat[(bestSplit[n]+1):n] <- bestAver[n]
   		n <- bestSplit[n]
 	}
 	return(yhat)
}


#Perform MAD winsorization:
madWins <- function(x,tau,k){
	xhat <- medianFilter(x,k)
	d <- x-xhat
	SD <- mad(d)
	z <- tau*SD
	xwin <- xhat + psi(d, z)
	outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
	return(list(ywin=xwin,sdev=SD,outliers=outliers))
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}





psi <- function(x,z){
 xwin <- x
 xwin[x < -z] <- -z
 xwin[x > z] <- z
 return(xwin)
}
