#Convert call-level to bout level and save

dir <- '~/Dropbox/vlad_callback_analysis/SunningCalls' #directory of files - only the files - for Ari
#dir <- 'C:/Users/Vlad/Dropbox/vlad_callback_analysis/data'
filename <- 'sunningoverlapanalysis_20180720.csv'
pb <- F #whether to use pb data instead

setwd(dir)

#LOAD DATA
if(!pb){
  data <- read.csv(filename,header=TRUE,stringsAsFactors=F) 
}

if(pb){
  data <- read.csv('/Users/astrandb/Dropbox/vlad_callback_analysis/data/LongformatAllplaybacks.csv',stringsAsFactors=F) #UNCOMMENT FOR PB
}

bout.level <- T

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

#Consistency for playbacks
if('Terminal' %in% colnames(data)){
  data$Terminal <- as.character(data$Terminal)
  data$terminal <- as.character(data$Terminal)
  data$terminal[which(data$terminal == 'Terminal')] <- 'terminal'
  data$Terminal[which(data$terminal == 'Terminal')] <- 'terminal'
}


#new way
data$prev.foc.call <- NA
max.inter.note.interval <- .25 #maximum time between notes to consider it still the same bout
if(bout.level){
  data$bout <- NA
  for(i in file.nums){
    curr.idxs <- which(data$file.id==i & data$label=='focal')
    bout.idx <- 1
    if(length(curr.idxs)==1){
      data$bout[curr.idxs[1]] <- bout.idx
    }
    if(length(curr.idxs)>1){
      for(j in seq(1,length(curr.idxs)-1,1)){
        data$bout[curr.idxs[j]] <- bout.idx
        gap <- data$start.time[curr.idxs[j+1]]-data$end.time[curr.idxs[j]]
        if(gap > max.inter.note.interval){
          bout.idx <- bout.idx + 1
        }
      }
    }
  }
}

#OLD WAY
#if(bout.level){
#  data$bout <- NA
#  for(i in file.nums){
#    curr.idxs <- which(data$file.id==i & data$label=='focal')
#    bout.idx <- 1
#    for(j in curr.idxs){
#      data$bout[j] <- bout.idx
#      if(data$terminal[j] == 'terminal'){
#        bout.idx <- bout.idx+1
#      }
#    }
#  }
#}

new.data <- data.frame()
for(i in 1:length(data.files.list)){
  curr <- data[which(data$FileID==data.files.list[i]),]
  curr$starttime <- curr$start.time
  curr$endtime <- curr$end.time
  curr.new <- data.frame()
  if(sum(curr$label=='focal')>0){
    bouts <- max(curr$bout,na.rm=T)
    if(bouts>0){
      for(j in 1:bouts){
        bout.idxs <- which(curr$bout==j)
        start <- min(curr$start.time[bout.idxs],na.rm=T)
        end <- max(curr$end.time[bout.idxs],na.rm=T)
        tmp <- curr[bout.idxs[1],]
        tmp$start.time <- start
        tmp$end.time <- end
        tmp$starttime <- start
        tmp$endtime <- end
        
        if(pb){
          startreal <- min(curr$Real.Start[bout.idxs],na.rm=T)
          endreal <- max(curr$RealEnd[bout.idxs],na.rm=T)
          tmp$Real.Start <- startreal
          tmp$RealEnd <- endreal
        }
        
        curr.new <- rbind(curr.new,tmp)
      }
    }
  }
  if(sum(curr$label!='focal')>0){
    curr.new <- rbind(curr.new,curr[which(curr$label!='focal'),])
  }
  curr.new <- curr.new[order(curr.new$start.time),]
  new.data <- rbind(new.data,curr.new)
}

data <- new.data

#clean up - any bouts longer than 10 sec are unrealistic (there was only 1 anyway)
data<- data[which((data$end.time-data$start.time)<10),]

if(!pb){
  save(list=c('data'),file='/Users/astrandb/Dropbox/vlad_callback_analysis/data/natural_sunning_data_bout_level.RData')
}

if(pb){
  save(list=c('data'),file='/Users/astrandb/Dropbox/vlad_callback_analysis/data/playback_data_bout_level.RData')
}

     