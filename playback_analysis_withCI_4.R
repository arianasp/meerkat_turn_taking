#Compute the focal call rate as a function of the width of window around a playback call
#Also look at inter-call and inter-bout interval distributions from control data

library(fitdistrplus)

use.natural.call.data <- F
randomize.times <- F
plot.CIs <- T
sample.by.file <- T

#LOAD DATA

if(use.natural.call.data){
	data <- read.csv('~/Dropbox/vlad_callback_analysis/SunningCalls/sunningoverlapanalysis_20180720.csv')
	data$filename <- data$FileID
	data$Real.Start <- data$starttime
	data$RealEnd <- data$endtime
} else{
	data <- read.table('~/Dropbox/vlad_callback_analysis/Data/LongformatAllplaybacks.csv',sep=',',header=T,stringsAsFactors=F)
}


files <- unique(data$filename)

windows.exp <- seq(-3,1.5,.1)
windows <- 10^(windows.exp) #time window (in sec)

call.rate.ctrl <- call.rate.pb <- rep(NA,length(windows))
call.rate.ctrl.all <- call.rate.pb.all <- files.all <- list()
for(k in 1:length(windows)){
	window <- windows[k]
	#loop over trials
	call.dur.tot.pb <- call.dur.tot.ctrl <- 0
	n.time.tot <- 0
	call.rate.ctrl.all[[k]] <- c(NA)
	call.rate.pb.all[[k]] <- c(NA)
	files.all[[k]] <- c(NA)
	for(i in 1:length(files)){
		file <- files[i]
		curr.dat <- data[which(data$filename==file),]
	
		if(use.natural.call.data){
			pb.dat <- curr.dat
			ctrl.dat <- curr.dat
		} else{
			#Separate playback and control data
			pb.dat <- curr.dat[which(curr.dat$Type=='Playback'),]
			ctrl.dat <- curr.dat[which(curr.dat$Type=='Control'),]
		}
		if((nrow(pb.dat)) > 0 & (nrow(ctrl.dat) > 0)){
			#Get times of playback calls
			if(use.natural.call.data){
				pb.call.idxs <- which(pb.dat$label=='background')
				pb.call.times <- pb.dat$RealEnd[pb.call.idxs] #calls are pinned at the end of each playback call
			} else{
				pb.call.idxs <- which(pb.dat$label=='playback')
				pb.call.times <- (pb.dat$Real.Start[pb.call.idxs] + pb.dat$RealEnd[pb.call.idxs])/2
				pb.call.times <- pb.dat$RealEnd[pb.call.idxs] #calls are pinned at the end of each playback call
			}
	
			#Separate out only calls from the focal
			pb.foc.dat <- pb.dat[which(pb.dat$label=='focal'),]
			ctrl.foc.dat <- ctrl.dat[which(ctrl.dat$label=='focal'),]
			
			if(use.natural.call.data){
				tmin.pb <- tmin.ctrl <- 0
				tmax.pb <- tmax.ctrl <- max(ctrl.dat$Real.Start,na.rm=T)
			} else{
				tmin.pb <- pb.dat$Real.Start[which(pb.dat$label=='Start Playback')]
				tmax.pb <- pb.dat$Real.Start[which(pb.dat$label=='End Playback')]
				tmin.ctrl <- ctrl.dat$Real.Start[which(ctrl.dat$label=='Start Control')]
				tmax.ctrl <- ctrl.dat$Real.Start[which(ctrl.dat$label=='End Control')]
			}
			
			if(randomize.times){
				pb.call.times <- runif(length(pb.call.times))*(tmax.pb - tmin.pb) + tmin.pb
			}
			
			
			#Deal with missing data
			if(length(tmin.pb)==0){
				tmin.pb <- 0
			}
			if(length(tmin.ctrl)==0){
				tmin.ctrl <- 0
			}
			if(length(tmax.pb)==0){
				tmax.pb <- max(pb.dat$RealEnd)
			}
			if(length(tmax.ctrl)==0){
				tmax.ctrl <- max(ctrl.dat$RealEnd)
			}
		
			tmin <- max(tmin.pb,tmin.ctrl)
			tmax <- min(tmax.pb,tmax.ctrl)
			#Get number of calls within a window of each playback call, as well as total amount of time in that window
			if(length(pb.call.times) > 0){
				for(j in 1:length(pb.call.times)){
					t0 <- pb.call.times[j] # - window/2
					tf <- pb.call.times[j] + window
					if((t0 > tmin) & (tf < tmax)){
						pb.call.idxs <- which((pb.foc.dat$Real.Start) <= tf & (pb.foc.dat$RealEnd >= t0))
						if(length(pb.call.idxs)>0){
							call.dur.tot.pb <- call.dur.tot.pb + sum(sapply(pb.foc.dat$RealEnd[pb.call.idxs],FUN=function(x){return(min(x,tf))})-sapply(pb.foc.dat$Real.Start[pb.call.idxs],FUN=function(x){return(max(x,t0))}))
						  call.rate.pb.all[[k]] <- c(call.rate.pb.all[[k]],sum(sapply(pb.foc.dat$RealEnd[pb.call.idxs],FUN=function(x){return(min(x,tf))})-sapply(pb.foc.dat$Real.Start[pb.call.idxs],FUN=function(x){return(max(x,t0))})) / (tf - t0))
						} else{
						  call.rate.pb.all[[k]] <- c(call.rate.pb.all[[k]],0)
						}
						ctrl.call.idxs <- which((ctrl.foc.dat$Real.Start) <= tf & (ctrl.foc.dat$RealEnd >= t0))
						if(length(ctrl.call.idxs)>0){
							call.dur.tot.ctrl <- call.dur.tot.ctrl + sum(sapply(ctrl.foc.dat$RealEnd[ctrl.call.idxs],FUN=function(x){return(min(x,tf))})-sapply(ctrl.foc.dat$Real.Start[ctrl.call.idxs],FUN=function(x){return(max(x,t0))}))
							call.rate.ctrl.all[[k]] <- c(call.rate.ctrl.all[[k]],sum(sapply(ctrl.foc.dat$RealEnd[ctrl.call.idxs],FUN=function(x){return(min(x,tf))})-sapply(ctrl.foc.dat$Real.Start[ctrl.call.idxs],FUN=function(x){return(max(x,t0))})) / (tf - t0))
						} else{
						  call.rate.ctrl.all[[k]] <- c(call.rate.ctrl.all[[k]],0)
						}
						files.all[[k]] <- c(files.all[[k]],file)
						n.time.tot <- n.time.tot + (tf - t0)
					}
				}
			}
		}
	}
	call.rate.ctrl[k] <- call.dur.tot.ctrl / n.time.tot * 100
	call.rate.pb[k] <- call.dur.tot.pb / n.time.tot * 100
}

#you may need to run it twice, once with randomize.times =T and once with randomize.times =F
if(randomize.times){
  call.rate.ctrl.rand <- call.rate.ctrl
  call.rate.ctrl.rand.all <- call.rate.ctrl.all
}

#bootstrap confidence intervals on means
if(plot.CIs){
  uppers.ctrl <- uppers.exp <- lowers.ctrl <- lowers.exp <- u75s.ctrl <- u75s.exp <- l25s.ctrl <- l25s.exp <- rep(NA,length(call.rate.ctrl.all))
  boots <- 1000
  for(i in 1:length(call.rate.ctrl.all)){
    if(use.natural.call.data){
      curr.exp <- call.rate.ctrl.all[[i]][2:length(call.rate.ctrl.all[[i]])] #start at 2 bc first is NA (from construction above)
      curr.ctrl <- call.rate.ctrl.rand.all[[i]][2:length(call.rate.ctrl.all[[i]])] #start at 2 bc first is NA (from construction above)
    } else{
      curr.exp <- call.rate.pb.all[[i]][2:length(call.rate.pb.all[[i]])] #start at 2 bc first is NA (from construction above)
      curr.ctrl <- call.rate.ctrl.all[[i]][2:length(call.rate.ctrl.all[[i]])] #start at 2 bc first is NA (from construction above)
    }
    curr.files <- files.all[[i]][2:length(files.all[[i]])] #start at 2 bc first is NA (from construction above)
    means.ctrl <- means.exp <- rep(NA,boots)
    for(j in 1:boots){
      
      #sample independently
      #exps <- sample(curr.exp,length(curr.exp),replace=T)
      #ctrls <- sample(curr.ctrl,length(curr.ctrl),replace=T)
      
      #block sampling (sample a whole trial)
      exps <- ctrls <- c()
      curr.files.unique <- unique(curr.files)
      for(k in 1:length(curr.files.unique)){
        f <- sample(curr.files.unique,1)
        exps <- c(exps,curr.exp[which(curr.files == f)])
        ctrls <- c(ctrls,curr.ctrl[which(curr.files == f)])
      }
      
      
      means.ctrl[j] <- mean(ctrls)*100
      means.exp[j] <- mean(exps)*100
    }
    uppers.ctrl[i] <- quantile(means.ctrl,0.975)
    uppers.exp[i] <- quantile(means.exp,0.975)
    lowers.ctrl[i] <- quantile(means.ctrl,0.025)
    lowers.exp[i] <- quantile(means.exp,0.025)
    u75s.ctrl[i] <- quantile(means.ctrl,0.75)
    u75s.exp[i] <- quantile(means.exp,0.75)
    l25s.ctrl[i] <- quantile(means.ctrl,0.25)
    l25s.exp[i] <- quantile(means.exp,0.25)
  }
}

if(use.natural.call.data){
  #quartz()
  quartz(width=6,height=6)
  par(mar=c(6,5,1,1))
	plot(NULL,xlim=c(min(windows),max(windows)),ylim=c(0,2.5),xlab='Time window (sec)',ylab='Call rate (% time calling)',log='x',cex.lab=2,cex.axis=2)
	
	if(plot.CIs){
	  polygon(c(windows,rev(windows)),c(uppers.ctrl,rev(lowers.ctrl)),col='#00000033',border=NA)
	  polygon(c(windows,rev(windows)),c(uppers.exp,rev(lowers.exp)),col='#00000066',border=NA)
	  #polygon(c(windows,rev(windows)),c(u75s.ctrl,rev(l25s.ctrl)),col='#00000066',border=NA)
	  #polygon(c(windows,rev(windows)),c(u75s.exp,rev(l25s.exp)),col='#FF000066',border=NA)
	}
	
	
	lines(windows,call.rate.ctrl.rand,col='black')
	lines(windows,call.rate.ctrl,lwd=1)
	legend(min(windows),2.5,c('Random times','After conspecific call'),cex=2,pch=c(17,19),lty=c(1,1))
	points(windows,call.rate.ctrl.rand,col='black',pch=17)
	points(windows,call.rate.ctrl,pch=19,bg='white')
} else{
  quartz(width=6,height=6)
  par(mar=c(6,5,1,1))
	plot(NULL,xlim=c(min(windows),max(windows)),ylim=c(0,2.5),xlab='Time window (sec)',ylab='Call rate (% time calling)',log='x',cex.lab=2,cex.axis=2)
	
	if(plot.CIs){
	  polygon(c(windows,rev(windows)),c(uppers.ctrl,rev(lowers.ctrl)),col='#00000033',border=NA)
	  polygon(c(windows,rev(windows)),c(uppers.exp,rev(lowers.exp)),col='#FF000066',border=NA)
	  #polygon(c(windows,rev(windows)),c(u75s.ctrl,rev(l25s.ctrl)),col='#00000066',border=NA)
	  #polygon(c(windows,rev(windows)),c(u75s.exp,rev(l25s.exp)),col='#FF000066',border=NA)
	}
	lines(windows,call.rate.ctrl,col='black')
	lines(windows,call.rate.pb,lwd=1)
	legend(min(windows),2.5,c('Control','Playback'),cex=2,pch=c(17,19),lty=c(1,1))
	points(windows,call.rate.ctrl,col='black',pch=17)
	points(windows,call.rate.pb,pch=19,bg='white')
}


#Statistical test - randomize labels. Test stat is difference in pct between ctrl and pb - not really valid because assumes independence where there isn't - don't use
n.rands <- 1000
test.stats.data <- unlist(lapply(call.rate.pb.all,function(x){return(mean(x>0,na.rm=T))},na.rm=T)) - unlist(lapply(call.rate.ctrl.all,function(x){return(mean(x>0,na.rm=T))},na.rm=T))


test.stats.rand <- matrix(NA,nrow=length(test.stats.data),ncol=n.rands)
for(r in 1:n.rands){
  for(i in 1:length(test.stats.data)){
    dat.curr <- c(call.rate.pb.all[[i]],call.rate.ctrl.all[[i]])
    #remove NAs
    dat.curr <- dat.curr[which(!is.na(dat.curr))]
    idxs.ctrl.rand <- sample(length(dat.curr),length(dat.curr)/2)
    idxs.pb.rand <- seq(1,length(dat.curr))[-idxs.ctrl.rand]
    pb.rand.curr <- dat.curr[idxs.pb.rand]
    ctrl.rand.curr <- dat.curr[idxs.ctrl.rand]
    test.stats.rand[i,r] <- mean(pb.rand.curr>0,na.rm=T) - mean(ctrl.rand.curr>0,na.rm=T)
  }
}

#get p values
ps <- rep(NA,length(test.stats.data))
for(i in 1:length(test.stats.data)){
  cumdist <- ecdf(test.stats.rand[i,])
  ps[i] <- cumdist(test.stats.data[i])
}
