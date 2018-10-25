#Code to determine the level of overlap in call sequences
#Does a time shift to determine whether this level changes when focal and background tracks are shifted relative to one another
#Also does a hypothesis test based on time shifts (not used in final paper)

#PARAMETERS

max.shift <- 0.5 #maximum value of the shifting window (in seconds)
n.rands <- 100 #number of randomizations to do
bout.level <- TRUE #whether to do analysis at the bout level (chunk bouts into all 1's) or call level
time.shift.randomization <- TRUE
synch.by.shift <- TRUE
other.plots <- FALSE


dir <- '~/Dropbox/vlad_callback_analysis/SunningCalls' #directory of files

setwd(dir)

#LOAD DATA

if(!bout.level){
  print('performing analysis at the note level')
  filename <- 'sunningoverlapanalysis_20180720.csv'
  print('reading in filename:')
  print(filename)
  data <- read.csv(filename,header=TRUE,stringsAsFactors=F)
}

if(bout.level){
  print('performing analysis at the bout level')
  print('reading in filename:')
  filename <- '~/Dropbox/vlad_callback_analysis/data/natural_sunning_data_bout_level.RData'
  load(filename)
}

#Rename columns appropriately if needed
if(is.null(data$FileID)){
	data$FileID <- data$filename
}

if(is.null(data$start.time)){
	data$start.time <- data$starttime
}

if(is.null(data$end.time)){
	data$end.time <- data$endtime
}

#Get file numbers
data.files.list <- unique(data$FileID)
file.nums <- 1:length(data.files.list)
data$file.id <- file.nums[match(data$FileID,data.files.list)]


#CONVERT TO BINARY STRINGS AND MAKE LISTS OF THESE STRINGS

focal.binary.strings <- background.binary.strings <- list()

for(j in 1:length(data.files.list)){
	curr <- data[which(data$file.id==j),] #get first file data
	seq.length <- max(curr$end.time)*1000 + max.shift*1000 #length of the vector of 1's and 0s indicating times of calls

	focal.calls <- rep(0,seq.length) #focal individual's calls
	background.calls <- rep(0,seq.length) #background calls

	for(i in 1:nrow(curr)){
		tstart <- curr$start.time[i]
		tend <- curr$end.time[i]
		
		if(bout.level){
			bout.id <- curr$bout[i]
			if(!is.na(bout.id)){
				tstart <- min(curr$start.time[which(curr$bout==bout.id)])
				tend <- max(curr$end.time[which(curr$bout==bout.id)])
			}
		}
		
		if(curr$label[i] %in% c('focal','Focal')){
			focal.calls[seq(tstart*1000,tend*1000)] <- 1
		}
		if(curr$label[i] %in% c('background','Background')){
			background.calls[seq(tstart*1000,tend*1000)] <- 1
		}
	
	}
	focal.binary.strings[[j]] <- focal.calls
	background.binary.strings[[j]] <- background.calls
	
}

#FUNCTIONS
#Compute the synch score for all pairs of strings
#Inputs:
#	focal.binary.strings: list of binary sequences for the focal individual (each element is a different file)
#	background.binary.strings: list of binary sequences for the background individuals (each element is a different file)
#Outputs:
#	tot.score: the synch score, computed across all sequences (between 0 and 1)
compute.synch.score <- function(focal.binary.strings,background.binary.strings){

	score <- tot.length <- 0
	for(i in 1:length(focal.binary.strings)){
		seq1 <- focal.binary.strings[[i]]
		seq2 <- background.binary.strings[[i]]
		score <- score + sum(seq1*seq2)
		tot.length <- tot.length + min(sum(seq1),sum(seq2))
	}

	tot.score <- score / tot.length 
	return(tot.score)
}

#Shift one of the sequences by an amount "shift"
#Inputs: 
#	seq1: binary vector of 1s and 0s
#	seq2: binary vector of 1s and 0s
#	max.shift: the maximum value to shift one of the seqs by (in seconds, defaults to 1)
#	shift: specifies how much to shift by (defaults to NULL, in which case the shift is chosen uniformly at random from the range 0 to max.shift)
# 	negative.shifts.allowed: whether negative shifts are allowed (in this case, sequences are randomly relabeled if a neg number is drawn)
#Outputs:
#	a matrix with the two sequences, one shifted, as columns
shift.seqs <- function(seq1,seq2,max.shift=1,shift=NULL,negative.shifts.allowed=T){
	#pick a random number to shift and shift sequence 1 by that number, remove the unmatched ends
	
	if(is.null(shift)){
		if(negative.shifts.allowed){
			shift <- round(runif(1,min=-max.shift,max=max.shift*1000))
		} else{
			shift <- round(runif(1,min=0,max=max.shift*1000))
		}
	}
	
	if(shift < 0){ #swap seqs
		curr <- seq1
		seq1 <- seq2
		seq2 <- curr
		shift <- -shift
	}
	seq1.cut <- seq1[seq(1,length(seq1)-shift)]
	seq1.shift <- rep(0,length(seq1))
	seq1.shift[seq(shift+1,length(seq1.cut)+shift)] <- seq1.cut
	seq1.shift <- seq1.shift[seq(shift+1,length(seq1.cut))]
	seq2.noshift <- seq2[seq(shift+1,length(seq1.cut))]
	return(cbind(seq1.shift,seq2.noshift))
}

#RANDOMIZATION TO COMPUTE SYNCH SCORES WITH RANDOM SHIFTS FROM 0 TO MAX.SHIFT
if(time.shift.randomization){
  print('performing time shift randomization')
  synch.scores.rand <- rep(NA,n.rands)
  for(i in 1:n.rands){
  	print(i)
  	focal.binary.strings.rand <- background.binary.strings.rand <- list()
  	for(j in 1:length(focal.binary.strings)){
  		new.seqs <- shift.seqs(focal.binary.strings[[j]],background.binary.strings[[j]],max.shift=max.shift)
  		focal.binary.strings.rand[[j]] <- new.seqs[,1]
  		background.binary.strings.rand[[j]] <- new.seqs[,2]
  	}
  	synch.scores.rand[i] <- compute.synch.score(focal.binary.strings.rand,background.binary.strings.rand)	
  }
  
  synch.score.data <- compute.synch.score(focal.binary.strings,background.binary.strings)
  
  #MAKE A HISTOGRAM PLOT
  histo <- hist(synch.scores.rand,breaks=seq(0,1,.001),plot=F)
  xmin <- histo$breaks[min(which(histo$counts !=0))-1]
  xmin <- min(xmin,synch.score.data-.005)
  xmax <- histo$breaks[max(which(histo$counts !=0))+1]
  
  plot(histo$mids,histo$counts/sum(histo$counts),type='l',xlim=c(xmin,xmax),pch=19,lwd=2,xlab='overlap rate',ylab='probability')
  points(histo$mids,histo$counts/sum(histo$counts),pch=19,cex=0.5)
  abline(v=synch.score.data,col='red',lwd=3,lty=2)
}

#Get synch score by time shift
if(synch.by.shift){
  print('computing synch score by time shift and plotting')
  shifts <- seq(-1000,1000,10)
  #shifts <- seq(-10,10,1)
  synch.scores.by.shift <- rep(NA,length(shifts))
  idx <- 1
  for(i in shifts){
  	print(i)
  	for(j in 1:length(focal.binary.strings)){
  		new.seqs <- shift.seqs(focal.binary.strings[[j]],background.binary.strings[[j]],max.shift=max.shift,shift=i)
  		focal.binary.strings.rand[[j]] <- new.seqs[,1]
  		background.binary.strings.rand[[j]] <- new.seqs[,2]
  	}
  	synch.scores.by.shift[idx] <- compute.synch.score(focal.binary.strings.rand,background.binary.strings.rand)
  	idx <- idx + 1
  }
  
  quartz(height=6,width=6)
  par(mar=c(6,6,1,1))
  plot(shifts/1000,synch.scores.by.shift,pch=19,cex=0.5,xlab='time shift (sec)',ylab='overlap rate',ylim=c(0,0.07),cex.axis=2,cex.lab=2,type='l',lwd=2)
  points(shifts/1000,synch.scores.by.shift,pch=19,cex=0.5)
  #lines(shifts/1000,synch.scores.by.shift,lwd=2)
  abline(v=0,col='red',lty=2,lwd=2)
}

#COMPUTE IBI AS A FUNCTION OF BCR
if(other.plots){
  print('making some other plots')
  window.size <- 10 #window size over which to compute BCR (in sec)
  
  n.bouts <- sum(!is.na(data$IBI)) #number of bouts
  bout.idxs <- which(!is.na(data$IBI)) #indexes to bouts in data frame
  
  bcrs <- ibis <- rep(NA,n.bouts) #vectors to hold IBIs and BCRs
  for(i in 1:length(bout.idxs)){
  	row <- bout.idxs[i]
  	ibis[i] <- data$IBI[row]
  	file <- data$FileID[row]
  	curr.dat <- data[which(data$FileID==file),]
  	curr.dat.bg <- curr.dat[which(curr.dat$label=='background'),]
  	if(nrow(curr.dat) == 0){
  		print(i)
  	}
  	tf <- data$end.time[row]
  	t0 <- tf - data$IBI[row]
  	mid <- (tf + t0) / 2
  	window.t0 <- mid - window.size/2
  	window.tf <- mid + window.size/2
  	if(window.t0 < 0 | window.tf > max(curr.dat$end.time,na.rm=T)){
  		bcrs[i] <- NA
  	} else{
  		calls <- sum(curr.dat.bg$start.time >= window.t0 & curr.dat.bg$end.time <= window.tf)
  		rate <- calls / window.size
  		bcrs[i] <- rate
  	}
  	if(nrow(curr.dat.bg)==0){
  		bcrs[i] <- NA
  	}
  }
  
  ibis.0 <- ibis[which(bcrs == 0)]
  ibis.non0 <- ibis[which(bcrs > 0)]
  
  bins <- c(0,quantile(bcrs,quants,na.rm=T),max(bcrs,na.rm=T))
  meds <- uppers <- lowers <- uppers2 <- lowers2 <- rep(NA,length(bins)-1)
  for(i in 1:(length(bins)-1)){
  	dat <- ibis[which(bcrs >= bins[i] & bcrs < bins[i+1])]
  	meds[i] <- median(dat,na.rm=T)
  	uppers[i] <- quantile(dat,0.75,na.rm=T)
  	lowers[i] <- quantile(dat,0.25,na.rm=T)
  	uppers2[i] <- quantile(dat,0.975,na.rm=T)
  	lowers2[i] <- quantile(dat,0.025,na.rm=T)
  }
  mids <- (bins[1:(length(bins)-1)] + bins[2:length(bins)])/2
  plot(NULL,xlim=c(0,max(mids)),ylim=c(0,max(uppers2,na.rm=T)),xlab='BCR (calls / sec)',ylab='IBI (sec)',main=paste('window size = ',window.size, ' sec',sep=''))
  polygon(c(mids,rev(mids)),c(uppers2,rev(lowers2)),col='#FF000022',border=NA)
  polygon(c(mids,rev(mids)),c(uppers,rev(lowers)),col='#FF000066',border=NA)
  lines(mids,meds,lwd=3)
  points(mids,meds,pch=19,cex=2)
  
  #Focal call rate vs. background call rate
  files <- unique(data$FileID)
  bcr <- fcr <- fbr <- rep(NA,length(files))
  for(i in 1:length(files)){
  	curr <- data[which(data$FileID==files[i]),]
  	t0 <- 0
  	tf <- max(data$end.time,na.rm=T)
  	dt <- tf - t0
  	bcr[i] <- sum(curr$label=='background',na.rm=T)/dt
  	fcr[i] <- sum(curr$label=='focal',na.rm=T)/dt
  	fbr[i] <- sum(curr$terminal=='terminal',na.rm=T)/dt
  }
  
  
  quants <- seq(0.2,0.8,0.2)
  bins <- c(0,quantile(bcr,quants,na.rm=T),max(bcr,na.rm=T))
  meds <- uppers <- lowers <- uppers2 <- lowers2 <- rep(NA,length(bins)-1)
  for(i in 1:(length(bins)-1)){
  	dat <- fbr[which(bcr >= bins[i] & bcr < bins[i+1])]
  	meds[i] <- median(dat,na.rm=T)
  	uppers[i] <- quantile(dat,0.75,na.rm=T)
  	lowers[i] <- quantile(dat,0.25,na.rm=T)
  	uppers2[i] <- quantile(dat,0.975,na.rm=T)
  	lowers2[i] <- quantile(dat,0.025,na.rm=T)
  }
  mids <- (bins[1:(length(bins)-1)] + bins[2:length(bins)])/2
  plot(NULL,xlim=c(min(mids),max(mids)),ylim=c(0,max(uppers2,na.rm=T)),xlab='Background Call Rate (calls / sec)',ylab='Focal Call Rate (calls / sec)')
  #points(bcr,fcr,pch=19,cex=0.3)
  polygon(c(mids,rev(mids)),c(uppers2,rev(lowers2)),col='#FF000022',border=NA)
  polygon(c(mids,rev(mids)),c(uppers,rev(lowers)),col='#FF000066',border=NA)
  lines(mids,meds,lwd=3)
  points(mids,meds,pch=19,cex=2)
  
  plot(NULL,xlim=c(min(bins),max(bins)),ylim=c(0,max(dat,na.rm=T)),xlab='Background Call Rate (calls / sec)',ylab='Focal Bout Rate (bouts / sec)')
  points(bcr,fbr,cex=0.5)
  for(i in 1:length(mids)){
  	polygon(c(bins[i],bins[i+1],bins[i+1],bins[i]),c(uppers2[i],uppers2[i],lowers2[i],lowers2[i]),col='#FF000022',border=NA)
  	polygon(c(bins[i],bins[i+1],bins[i+1],bins[i]),c(uppers[i],uppers[i],lowers[i],lowers[i]),col='#FF000066',border=NA)
  	lines(c(bins[i],bins[i+1]),c(meds[i],meds[i]),lwd=2)
  }
}
  








