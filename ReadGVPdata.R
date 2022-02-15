a <- dir('DataFolder')
ndatset <- length(a)

###############################################
## For all the data available in GVP records
for (kk in 1:ndatset)
{
data <- read.csv(paste("DataFolder/",a[kk],sep=""),header=T)
data1 <- data
data$Start.Month[is.na(data$Start.Month)] <- 0
data$Start.Day[is.na(data$Start.Day)] <- 0
# use July and 15th for the missing months and days for both the 
# historical and non-historical records. 

useM <- data$Start.Month==0
useD <- data$Start.Day==0
indM <- (1:nrow(data))[useM]
indD <- (1:nrow(data))[useD]
for (iM in indM){
	data$Start.Month[iM] <- 7
}
mth <- data$Start.Month
for (iD in indD){
  data$Start.Day[iD] <- 15
}

# Convert year.month.day to days
t <- NULL
sty <- data$Start.Year 
mth <- data$Start.Month
std <- data$Start.Day
for (i in 1:length(sty)){
  if (sty[i] >= 0){
    t[i] <- sty[i]+
         (julian(as.Date(paste(sty[i],mth[i],std[i],sep="-")))-
           julian(as.Date(paste(sty[i]-1,12,31,sep="-"))))/
           (julian(as.Date(paste(sty[i],12,31,sep="-")))-
           julian(as.Date(paste(sty[i]-1,12,31,sep="-"))))
  }else{ # make all negative years positive so julian will work
    t[i] <- sty[i]+
         (julian(as.Date(paste(sty[i]+10000,mth[i],std[i],sep="-")))-
           julian(as.Date(paste(sty[i]-1+10000,12,31,sep="-"))))/
           (julian(as.Date(paste(sty[i]+10000,12,31,sep="-")))-
           julian(as.Date(paste(sty[i]-1+10000,12,31,sep="-"))))
  }
}

# For duplicate dates, the same year but with no month or day given 
# for the duplicates, divide the year evenly among the duplicate dates
tt <- t
ind <- 1:length(tt)
for (j in 1:length(tt)){
  ndup <- sum(t[-j]==t[j]) # number of duplicates
  if (ndup>=0.5){
    temind <- ind[t[-j]==t[j]]
    tt[c(j,temind)] <- t[j]+1:(ndup+1)/(ndup+2)
  } 
}

dat <- data
dat$ti <- tt
datt <- dat[order(dat$ti),]
eval(parse(text=paste('d',substr(a[kk],1,4),'=datt',sep="")))
eval(parse(text=paste('ti',substr(a[kk],1,4),'=d',substr(a[kk],1,4),'$ti',sep="")))

###############################################
## For historical observations only
datth <- datt[datt$Evidence.Method..dating.=="Historical Observations",]
eval(parse(text=paste('hisd',substr(a[kk],1,4),'=datth',sep="")))
eval(parse(text=paste('histi',substr(a[kk],1,4),'=hisd',substr(a[kk],1,4),'$ti',sep="")))

###############################################
## For VEI 3 only
datt$VEI[is.na(datt$VEI)] <- 4.5
datt3 <- datt[datt$VEI>=2.99,]
eval(parse(text=paste('dV3',substr(a[kk],1,4),'=datt3',sep="")))
eval(parse(text=paste('tiV3',substr(a[kk],1,4),'=dV3',substr(a[kk],1,4),'$ti',sep="")))

###############################################
## For VEI 4 only
datt$VEI[is.na(datt$VEI)] <- 4.5
datt4 <- datt[datt$VEI>=3.99,]
eval(parse(text=paste('dV4',substr(a[kk],1,4),'=datt4',sep="")))
eval(parse(text=paste('tiV4',substr(a[kk],1,4),'=dV4',substr(a[kk],1,4),'$ti',sep="")))
}


#get(paste('d',substr(a[kk],1,4),sep=""))
#get(paste('ti',substr(a[kk],1,4),sep=""))

#get(paste('hisd',substr(a[kk],1,4),sep=""))
#get(paste('histi',substr(a[kk],1,4),sep=""))


#cbind(data1$Start.Year,data1$Start.Month,data1$Start.Day,dat$ti)








