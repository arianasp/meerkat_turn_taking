#Background call rate vs focal call rate

filename <- '/Users/astrandb/Dropbox/vlad_callback_analysis/SunningCalls/sunningoverlapanalysis_20180720.csv'

data <- read.csv(filename)

files <- unique(data$FileID)
foc.calls <- bg.calls <- times <- rep(NA,length(files))

for(i in 1:length(files)){
  
  curr <- data[which(data$FileID==files[i]),]
  max.t <- max(curr$endtime)
  
  n.foc.calls <- sum(curr$label=='focal')
  n.bg.calls <- sum(curr$label=='background')
  
  times[i] <- max.t
  foc.calls[i] <- n.foc.calls
  bg.calls[i] <- n.bg.calls
  
  
}

bg.call.rates <- bg.calls/times*60
foc.calls.rates <-  foc.calls/times*60

bins <- quantile(bg.call.rates,seq(0,1,1/6))

medians <- uppers <- lowers <- uppers2 <- lowers2 <- rep(NA,length(bins)-1)

for(i in 1:(length(bins)-1)){
  curr.idxs <- which((bg.call.rates >= bins[i]) & (bg.call.rates < bins[i+1]))
  medians[i] <- median(foc.calls.rates[curr.idxs])
  uppers[i] <- quantile(foc.calls.rates[curr.idxs],0.75)
  lowers[i] <- quantile(foc.calls.rates[curr.idxs],0.25)
  uppers2[i] <- quantile(foc.calls.rates[curr.idxs],0.975)
  lowers2[i] <- quantile(foc.calls.rates[curr.idxs],0.025)
}

quartz(width=6,height=6)
par(mar=c(6,5,1,1))
plot(NULL,xlab='Background call rate (calls / min)',ylab='Focal call rate (calls / min)',cex.lab=2,cex.axis=2,pch=19,xlim=c(0,max(bg.call.rates)),ylim=c(0,max(foc.calls.rates)))
for(i in 1:length(medians)){
  polygon(c(bins[i],bins[i+1],bins[i+1],bins[i]),c(uppers2[i],uppers2[i],lowers2[i],lowers2[i]),col='#FF000033',border=NA)
  polygon(c(bins[i],bins[i+1],bins[i+1],bins[i]),c(uppers[i],uppers[i],lowers[i],lowers[i]),col='#FF000066',border=NA)
  lines(c(bins[i],bins[i+1]),c(medians[i],medians[i]),lwd=4)
}
points(bg.call.rates,foc.calls.rates,cex=1,lwd=2)

