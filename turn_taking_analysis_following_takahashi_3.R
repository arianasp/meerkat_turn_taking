#Turn-taking analysis following marmoset paper
#Takahashi 2013, "Coupled oscillator dynamics of vocal turn-taking in monkeys"


#-------------------NATURAL CALLS---------------
#Load data
#load('/Users/astrandb/Dropbox/vlad_callback_analysis/data/natural_sunning_data_bout_level.RData') #Uncomment for bout level analysis
data <- read.csv('/Users/astrandb/Dropbox/vlad_callback_analysis/SunningCalls/sunningoverlapanalysis_20180720.csv') #Uncomment for individual level analysis

#Autocorrelation (distribution of inter-bout intervals for focal) - compare to Figure 3D in Marmoset paper
#Upshot: Eyeballing it, the calls appear not to be periodic at the individual level (don't see bumps in prob distrib at integer intervals)
foc.dat <- data[which(data$label=='focal'),]
ibis <- c()
files <- unique(foc.dat$FileID)
for(i in 1:length(files)){
  curr <- foc.dat[which(foc.dat$FileID==files[i]),]
  if(nrow(curr)>=2){
    for(j in 1:nrow(curr)){
      ibis <- c(ibis,curr$starttime[(j+1):nrow(curr)]-curr$endtime[j])
    }
    #ibis <- c(ibis,curr$starttime[2:nrow(curr)] - curr$endtime[1:(nrow(curr)-1)])
  }
}

ibis <- ibis[!is.na(ibis)]
dens <- density(ibis,n=1000,from=0,to=600)
quartz(width=6,height=6)
par(mar=c(6,5,1,1))
plot(dens$x,dens$y,type='l',lwd=5,col='red',xlab='Inter-bout interval (sec)',ylab='Probability density',xlim=c(0,300),cex.lab=2,cex.axis=2)

#Get the distribution of bout lengths (for focal) and call lengths (for background)
dens.foc <- density(data$duration[which(data$label=='focal')],n=1000,from=0,to=5)
dens.bg <- density(data$duration[which(data$label=='background')],n=1000,from=0,to=5)
plot(dens.foc$x,dens.foc$y,type='l')
lines(dens.bg$x,dens.bg$y,col='red')

#FIGURE 2 FROM MARMOSET PAPER
#Check 'reset' hypothesis (compare to Figure 2 of marmoset paper)
#Compute inter-bout interval when there was an intervening background call vs not
ibi.bf <- ibi.ff <- c()
files <- unique(data$FileID)
data$label <- as.character(data$label) 
for(i in 1:length(files)){
  #get current data
  curr <- data[which(data$FileID==files[i]),]
  
  #make sure sorted by time
  curr <- curr[order(curr$starttime),]
  
  if(nrow(curr)>=3){
    curr$prev.lab <- curr$prev.end <- curr$prevprev.lab <- NA
    curr$prev.lab[2:nrow(curr)] <- curr$label[1:(nrow(curr)-1)]
    curr$prev.end[2:nrow(curr)] <- curr$endtime[1:(nrow(curr)-1)]
    curr$prevprev.lab[3:nrow(curr)] <- curr$label[1:(nrow(curr)-2)]
    curr$ibi <-  curr$starttime - curr$prev.end

    #Variant 1: Use all intervals of type background-focal and all of type focal-focal
    ibi.bf <- c(ibi.bf,curr$ibi[which(curr$label=='focal' & curr$prev.lab=='background')])
    ibi.ff <- c(ibi.ff,curr$ibi[which(curr$label=='focal' & curr$prev.lab=='focal')])
    
    #Variant 2: Use all intervals of type focal-focal BUT only intervals of type background-focal for which the preceding call was focal
    #(i.e. only one background call heard between 2 focal calls)
    #ibi.bf <- c(ibi.bf,curr$ibi[which(curr$label=='focal' & curr$prev.lab=='background' & curr$prevprev.lab=='focal')])
    #ibi.ff <- c(ibi.ff,curr$ibi[which(curr$label=='focal' & curr$prev.lab=='focal')])

    }
}

#Check 'inhibition' hypothesis - swap background and focal from different files, remove overlaps, then look at inter-bout intervals of focal
#Also compute overlap rate and overlaps per minute in randomized data
n.rands <- 100
pvals <- rep(NA,n.rands)
ovlps.rand <- ts.rand <- ovlp.time.rand <- max.ovlp.time.rand <- rep(NA,n.rands)
ibi.inhib.all <- ibi.ff.all <- list()
#options(warn=2)
for(pp in 1:n.rands){
  print(pp)
  #First generate permuted dataset (with overlaps removed) and compute IBIs for focal there
  files <- unique(data$FileID)
  files.perm <- sample(files)
  #files.perm <- files
  data.foc <- data[which(data$label=='focal'),]
  data.bg <- data[which(data$label=='background'),]
  data.perm <- data.frame()
  ovlps.rand[pp] <- 0
  ts.rand[pp] <- 0
  ovlp.time.rand[pp] <- 0
  max.ovlp.time.rand[pp] <- 0
  for(i in 1:length(files)){
    curr.foc <- data.foc[which(data.foc$FileID==files[i]),]
    curr.bg <- data.bg[which(data.bg$FileID==files.perm[i]),]
    other.foc <- data.foc[which(data.foc$FileID==files.perm[i]),]
    other.bg <- data.bg[which(data.bg$FileID==files[i]),]
    max.t <- min(max(c(curr.foc$endtime,other.bg$endtime)),max(c(curr.bg$endtime,other.foc$endtime)))
    
    curr.foc <- curr.foc[which(curr.foc$endtime <= max.t),]
    curr.bg <- curr.bg[which(curr.bg$endtime <= max.t),]
    
    if(nrow(curr.foc)>0 & nrow(curr.bg)>0){
      max.ovlp.time.rand[pp] <- max.ovlp.time.rand[pp] + min(sum(curr.foc$endtime-curr.foc$starttime),sum(curr.bg$endtime-curr.bg$starttime))
      curr.bg$FileID <- rep(files[i],nrow(curr.bg))
      ts.rand[pp] <- ts.rand[pp] + max.t
      curr.foc.ovlp.rem <- data.frame()
      for(j in 1:nrow(curr.foc)){
        ovlps <- which((curr.bg$starttime <= curr.foc$endtime[j]) & (curr.bg$endtime >= curr.foc$starttime[j]))
        if(length(ovlps)>0){
          ovlps.rand[pp] <- ovlps.rand[pp]+1
          for(o in ovlps){
            ovlp.time.rand[pp] <- ovlp.time.rand[pp] + min(curr.bg$endtime[o],curr.foc$endtime[j]) - max(curr.bg$starttime[o],curr.foc$starttime[j])
          }
          #if the focal call starts first, remove background calls, otherwise remove (actually, don't add) focal calls
          if(sum(curr.bg$starttime[ovlps] < curr.foc$starttime[j])==0){
            curr.foc.ovlp.rem <- rbind(curr.foc.ovlp.rem,curr.foc[j,])
            curr.bg <- curr.bg[-ovlps,]
          }
        } else{ #if no overlap, include focal call
          curr.foc.ovlp.rem <- rbind(curr.foc.ovlp.rem,curr.foc[j,])
        }
      }
      data.perm <- rbind(data.perm,curr.foc.ovlp.rem)
      data.perm <- rbind(data.perm,curr.bg)
    }
  }
  
  #Then get intervals between focal-focal sequences
  ibi.inhib <- c()
  for(i in 1:length(files)){
    curr <- data.perm[which(data.perm$FileID==files[i]),]
    curr <- curr[order(curr$starttime),]
    if(nrow(curr)>=2){
      curr$prev[2:nrow(curr)] <- curr$label[1:(nrow(curr)-1)]
      curr$ibi[2:nrow(curr)] <- curr$starttime[2:nrow(curr)] - curr$endtime[1:(nrow(curr)-1)]
      ibi.inhib <- c(ibi.inhib,curr$ibi[which(curr$label=='focal' & curr$prev=='focal')])
    }
  }
  
  ibi.inhib.all[[pp]] <- ibi.inhib
  ibi.ff.all[[pp]] <- ibi.ff
  
  #pvals[pp] <- ks.test(ibi.inhib,ibi.ff,alternative='greater')$p
}

#Get overlap rate in the real data
ovlps.data <- ts.data <- ovlp.time.data <- max.ovlp.time.data <- 0
options(warn=1)
files <- unique(data$FileID)
data.foc <- data[which(data$label=='focal'),]
data.bg <- data[which(data$label=='background'),]
for(i in 1:length(files)){
  curr.foc <- data.foc[which(data.foc$FileID==files[i]),]
  curr.bg <- data.bg[which(data.bg$FileID==files[i]),]
  max.t <- max(max(curr.foc$endtime),max(curr.bg$endtime))
    
  curr.foc <- curr.foc[which(curr.foc$endtime <= max.t),]
  curr.bg <- curr.bg[which(curr.bg$endtime <= max.t),]
  
  if(nrow(curr.foc)>0 & nrow(curr.bg)>0){
    max.ovlp.time.data <- max.ovlp.time.data + min(sum(curr.foc$endtime-curr.foc$starttime),sum(curr.bg$endtime-curr.bg$starttime))
    ts.data <- ts.data + max.t
    for(j in 1:nrow(curr.foc)){
      ovlps <- which((curr.bg$starttime <= curr.foc$endtime[j]) & (curr.bg$endtime >= curr.foc$starttime[j]))
      if(length(ovlps)>0){
        ovlps.data <- ovlps.data+1
        for(o in ovlps){
          ovlp.time.data <- ovlp.time.data + min(curr.bg$endtime[o],curr.foc$endtime[j]) - max(curr.bg$starttime[o],curr.foc$starttime[j])
        }
      }
    }
  }
}


#Do overlaps occur less often than in randomized dataset?
quartz(height=6,width=6)
par(mar=c(6,5,1,1))
hist(ovlp.time.rand/max.ovlp.time.rand,breaks=seq(0,1,.001),xlab='Overlap rate',ylab='Frequency',cex.lab=2,cex.axis=2,main='',col='gray',xlim=c(0,.08)) #overlap rate histogram (ovlps per min)
abline(v=ovlp.time.data/max.ovlp.time.data,lty=2,col='red',lwd=3)

#Get density distributions
n.samps <- 10000
dens.bf <- density(ibi.bf,n=n.samps,from=0,to=60)
dens.ff <- density(ibi.ff,n=n.samps,from=0,to=60)
dens.inhib <- density(ibi.inhib,n=n.samps,from=0,to=60)


#Same thing but with bootstrapped confidence intervals
n.boots <- 100
bf.y <- ff.y <- inhib.y <- matrix(NA,nrow=n.boots,ncol=length(dens.bf$x))
for(r in 1:n.boots){
  print(r)
  dens.bf.boot <- density(sample(ibi.bf,replace=T),n=n.samps,from=0,to=60)
  dens.ff.boot <- density(sample(ibi.ff,replace=T),n=n.samps,from=0,to=60)
  #dens.inhib.boot <- density(sample(ibi.inhib,replace=T),n=n.samps,from=0,to=60)
  dens.inhib.boot <- density(sample(ibi.inhib.all[[r]],replace=T),n=n.samps,from=0,to=60)
  bf.y[r,] <- dens.bf.boot$y
  ff.y[r,] <- dens.ff.boot$y
  inhib.y[r,] <- dens.inhib.boot$y
}

bf.uppers <- apply(bf.y,2,function(x){return(quantile(x,0.975))})
bf.lowers <- apply(bf.y,2,function(x){return(quantile(x,0.025))})
bf.medians <- apply(bf.y,2,function(x){return(quantile(x,.5))})
ff.uppers <- apply(ff.y,2,function(x){return(quantile(x,0.975))})
ff.lowers <- apply(ff.y,2,function(x){return(quantile(x,0.025))})
ff.medians <- apply(ff.y,2,function(x){return(quantile(x,0.5))})
inhib.uppers <- apply(inhib.y,2,function(x){return(quantile(x,0.975))})
inhib.lowers <- apply(inhib.y,2,function(x){return(quantile(x,0.025))})
inhib.medians <- apply(inhib.y,2,function(x){return(quantile(x,0.5))})
ymax <- max(cbind(bf.uppers,ff.uppers,inhib.uppers))

#Make plot
idxs1 <- which(dens.bf$x>=0)
idxs2 <- which(dens.ff$x>=0)
idxs3 <- which(dens.inhib$x>=0)
col1 <- '#008800'
col2 <- '#888888'
col3 <- '#0000FF'
quartz(width=6,height=6)
par(mar=c(6,5,1,1))
plot(NULL,xlim=c(0,1),ylim=c(0.01,ymax),xlab='Interval (sec)',ylab='Probability density',cex.lab=2,cex.axis=2)
polygon(c(dens.ff$x[idxs2],rev(dens.ff$x[idxs2])),c(ff.uppers,rev(ff.lowers)),col=paste(col1,'66',sep=''),border=NA)
polygon(c(dens.bf$x[idxs1],rev(dens.bf$x[idxs1])),c(bf.uppers,rev(bf.lowers)),col=paste(col2,'66',sep=''),border=NA)
polygon(c(dens.inhib$x[idxs3],rev(dens.inhib$x[idxs3])),c(inhib.uppers,rev(inhib.lowers)),col=paste(col3,'66',sep=''),border=NA)

lines(dens.ff$x[idxs2],ff.medians[idxs2],type='l',col=col1,lwd=4)
lines(dens.bf$x[idxs1],bf.medians[idxs1],type='l',col=col2,lwd=4)
lines(dens.inhib$x[idxs3],inhib.medians[idxs3],type='l',col=col3,lwd=4)
legend('topright',lty=c(1,1,1),col=c(col1,col2,col3),legend=c('Data','Reset','Inhibition'),lwd=c(5,5),cex=1.5)

#stats
ks.test(ibi.bf,ibi.ff) #D = 0.318, p < 0.001 (bout level)
ks.test(ibi.inhib,ibi.ff) #D = 0.051 p = 0.63 (bout level) / D = 0.095 p = 0.03 (variant 2) why is ff different?

#------------------PLAYBACK EXPERIMENT----------------
n.samps <- 10000
data <- read.csv('/Users/astrandb/Dropbox/vlad_callback_analysis/data/LongformatAllplaybacks.csv') #call level
#load('/Users/astrandb/Dropbox/vlad_callback_analysis/data/playback_data_bout_level.RData')
data.bouts <- data
#get inter-bout intervals for playback - focal (in PB) and focal - focal (in both)

#Get intervals between playback and focal, and between focal and focal
ibi.ff.pb <- ibi.pf <- c()
files <- unique(data.bouts$filename)
for(i in 1:length(files)){
  curr <- data.bouts[which(data.bouts$filename==files[i]),]
  curr$prev.lab <- NA
  curr$prev.lab[2:nrow(curr)] <- as.character(curr$label[1:(nrow(curr)-1)])
  curr$prev.t[2:nrow(curr)] <- curr$start.time[2:nrow(curr)] - curr$end.time[1:(nrow(curr)-1)]
  ibi.pf <- c(ibi.pf,curr$prev.t[which(curr$label=='focal' & curr$prev.lab=='playback' & curr$Type=='Playback')])
  ibi.ff.pb <- c(ibi.ff.pb,curr$prev.t[which(curr$label=='focal' & curr$prev.lab=='focal' & curr$Type=='Playback')])
}

#calculate density
dens.ff.pb <- density(ibi.ff.pb,n=n.samps,from=0,to=60)
dens.pf <- density(ibi.pf,n=n.samps,from=0,to=60)


#Randomization - swap control for playback, remove overlaps, then recompute F-F intervals
data.perm <- data.frame()
for(i in 1:length(files)){
  curr <- data.bouts[which(data.bouts$filename==files[i]),]
  curr.pb <- curr[which(curr$Type=='Playback' & curr$label=='playback'),]
  curr.ctrl <- curr[which(curr$Type=='Control'),]
  curr.ctrl <- rbind(curr.ctrl,curr.pb)
  curr.ctrl <- curr.ctrl[order(curr.ctrl$Real.Start),]
  tmax <- min(max(curr.ctrl$RealEnd[which(curr.ctrl$label=='focal')]),max(curr.ctrl$RealEnd[which(curr.ctrl$label=='playback')]))
  dat.swap <- curr.ctrl[which(curr.ctrl$RealEnd <= tmax),]
  
  #find and get rid of overlapping intervals
  curr.foc <- dat.swap[which(dat.swap$label=='focal'),]
  curr.bg <- dat.swap[which(dat.swap$label=='playback'),]
  curr.foc.ovlp.rem <- data.frame()
  if(nrow(curr.foc)>0 & nrow(curr.bg)>0){
    for(j in 1:nrow(curr.foc)){
      ovlps <- which((curr.bg$Real.Start <= curr.foc$RealEnd[j]) & (curr.bg$RealEnd >= curr.foc$Real.Start[j]))
      if(length(ovlps)>0){
        #if the focal call starts first, remove background calls, otherwise remove (actually, don't add) focal calls
        if(sum(curr.bg$Real.Start[ovlps] < curr.foc$Real.Start[j])==0){
          curr.foc.ovlp.rem <- rbind(curr.foc.ovlp.rem,curr.foc[j,])
          curr.bg <- curr.bg[-ovlps,]
        }
      } else{ #if no overlap, include focal call
        curr.foc.ovlp.rem <- rbind(curr.foc.ovlp.rem,curr.foc[j,])
      }
    }
    data.perm <- rbind(data.perm,curr.foc.ovlp.rem)
    data.perm <- rbind(data.perm,curr.bg)
  }
}

#Then get intervals between focal-focal sequences
ibi.inhib <- c()
for(i in 1:length(files)){
  curr <- data.perm[which(data.perm$filename==files[i]),]
  curr <- curr[order(curr$Real.Start),]
  if(nrow(curr)>=2){
    curr$prev[2:nrow(curr)] <- as.character(curr$label[1:(nrow(curr)-1)])
    curr$ibi[2:nrow(curr)] <- curr$Real.Start[2:nrow(curr)] - curr$RealEnd[1:(nrow(curr)-1)]
    print(length(which(curr$label=='focal' & curr$prev=='focal')))
    ibi.inhib <- c(ibi.inhib,curr$ibi[which(curr$label=='focal' & curr$prev=='focal')])
  }
}

dens.inhib <- density(ibi.inhib,n=n.samps,from=0,to=60)
ymax <-max(max(dens.ff.pb$y,dens.inhib$y),dens.pf$y)


#Confidence intervals bootstrapped
n.boots <- 1000
pf.y <- ff.y <- inhib.y <- matrix(NA,nrow=n.boots,ncol=length(dens.pf$x))
for(r in 1:n.boots){
  print(r)
  dens.pf.boot <- density(sample(ibi.pf,replace=T),n=n.samps,from=0,to=60)
  dens.ff.boot <- density(sample(ibi.ff.pb,replace=T),n=n.samps,from=0,to=60)
  dens.inhib.boot <- density(sample(ibi.inhib,replace=T),n=n.samps,from=0,to=60)
  pf.y[r,] <- dens.pf.boot$y
  ff.y[r,] <- dens.ff.boot$y
  inhib.y[r,] <- dens.inhib.boot$y
}

pf.uppers <- apply(pf.y,2,function(x){return(quantile(x,0.975))})
pf.lowers <- apply(pf.y,2,function(x){return(quantile(x,0.025))})
ff.pb.uppers <- apply(ff.y,2,function(x){return(quantile(x,0.975))})
ff.pb.lowers <- apply(ff.y,2,function(x){return(quantile(x,0.025))})
inhib.uppers <- apply(inhib.y,2,function(x){return(quantile(x,0.975))})
inhib.lowers <- apply(inhib.y,2,function(x){return(quantile(x,0.025))})
ymax <- max(cbind(pf.uppers,ff.pb.uppers,inhib.uppers))

#colors
col1 <- '#008800'
col2 <- '#888888'
col3 <- '#0000FF'

#plot
quartz(width=6,height=6)
par(mar=c(6,5,1,1))
plot(NULL,xlim=c(0,2),ylim=c(0,ymax),xlab='Interval (sec)',ylab='Probability density',cex.lab=2,cex.axis=2)
polygon(c(dens.ff.pb$x,rev(dens.ff.pb$x)),c(ff.pb.uppers,rev(ff.pb.lowers)),col=paste(col1,'66',sep=''),border=NA)
polygon(c(dens.pf$x,rev(dens.pf$x)),c(pf.uppers,rev(pf.lowers)),col=paste(col2,'66',sep=''),border=NA)
polygon(c(dens.inhib$x,rev(dens.inhib$x)),c(inhib.uppers,rev(inhib.lowers)),col=paste(col3,'66',sep=''),border=NA)

lines(dens.ff.pb$x,dens.ff.pb$y,type='l',col=col1,lwd=4)
lines(dens.pf$x,dens.pf$y,type='l',col=col2,lwd=4)
lines(dens.inhib$x,dens.inhib$y,type='l',col=col3,lwd=4)

legend('topright',lty=c(1,1,1),lwd=c(5,5,5),col=c(col1,col2,col3),legend=c('Data','Reset','Inhibition'),cex=1.5)

#stats
ks.test(ibi.inhib,ibi.ff.pb) #D = 0.118, p = 0.93 (2 sided)
ks.test(ibi.pf,ibi.ff.pb) #D = 0.365, p = 0.001 (2 sided)


#WHEN TOGETHER, ARE INTERVALS BETWEEN FOCAL BOUTS SHOWING PERIODICITY?
ibi.foc <- c()
for(i in 1:length(files)){
  print(i)
  curr <- data[which(data$FileID==files[i]),]
  if((nrow(curr)>2)){
    for(j in 1:(nrow(curr)-1)){
      ibi.foc <- c(ibi.foc,curr$starttime[(j+1):nrow(curr)] - curr$endtime[j])
    }
  }
}

n.samps <- 10000
dens.ibi.foc <- density(ibi.foc,n=n.samps,from=0,to=1200)

quartz(width=6,height=6)
par(mar=c(6,5,1,1))
plot(dens.ibi.foc$x,dens.ibi.foc$y,xlim=c(0,400),type='l',lwd=4,col='red',xlab='Interval (sec)',ylab='Probability density',cex.lab=2,cex.axis=2,ylim=c(0,.01))

