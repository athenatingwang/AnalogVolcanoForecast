## VEI 3+
source("ReadGVPdata.R")

## Using eruption records with minimum VEI 3

aname <- c(
"Tongariro",
"Asamayama",             
"Avachinsky",           
"Azumayama",                            
"Bezymianny",
"Calbuco",              
"Chokaisan",             
"Colima",               
"Fujisan",               
"Gamalama",             
"Guagua Pichincha",      
"Kuchinoerabujima",                  
"Miyakejima",            
"Nevado del Ruiz",        
"Pico de Orizaba",         
"Popocatepetl",         
"Puyehue-Cordon Caulle",  
"Sheveluch",            
"Suwanosejima",          
"Tokachidake",          
"Tolbachik",             
"Tungurahua",            
"Villarrica")  

Nind1 <- c(21,1:20,22:23)
Nind <- Nind1[1]

## Catalog last accessed on 18 September 2018
## Convert the last access date to days in 2018
lastday <- 2018+(julian(as.Date(paste(2018,9,18,sep="-")))-
           julian(as.Date(paste(2017,12,31,sep="-"))))/
           (julian(as.Date(paste(2018,12,31,sep="-")))-
           julian(as.Date(paste(2017,12,31,sep="-"))))

##### Obtain inter-event times for each volcano, paste(substr(a[kj],1,4),'IE',sep="")
## temp is the times of volcanic eruptions
## max.erp is the time of the last eruption in the GVP catalog
##### Obtain the censored time for each volcano, paste("y",substr(a[kj],1,4)
## len.volc is the number of events from each volcano

max.erp <- len.volc <- NULL
for (kj in 1:ndatset){
  ## tiV3: VEI 3+ eriptions only
  temp <- eval(parse(text=paste("tiV3",substr(a[kj],1,4),sep="")))
  len.volc[kj] <- length(temp)
  max.erp <- c(max.erp,temp[length(temp)])
  eval(parse(text=paste(substr(a[kj],1,4),'IE=diff(temp)/1000',sep="")))
  eval(parse(text=paste("y",substr(a[kj],1,4),'=(lastday[1]-temp[length(temp)])/1000',sep="")))
}

## Combine all data together, making sure that Tongariro is the first
## inter.censor contains the inter-event times plus last censored time as NA
##   from all the volcanoes, with Tongariro as the first
## volclen is the number of events from each volcano, with Tongariro as the first
## repi is a vector that contains the repeated index of each volcano for
##   the length of the number of eruptions from each volcano
## rep01 is a vector that contains repeated 0's for interevent times and 1 for
##   censored times
## The fourth column of dat (Interval) contains the 
##   maximum interevent times+1 for 
##   non-censored data and the censored time from the last censored time

inter.censor <- volclen <- repi <- rep01 <- Interval <- NULL
N2tem <- 0 # number of events from each volcano
i <- 0
for (jk in Nind){
  tem1 <- get(paste(substr(a[jk],1,4),'IE',sep=""))
  tem2 <- get(paste("y",substr(a[jk],1,4),sep=""))
  tem <- c(tem1,NA)
  max.inter <- max(tem1)+1
  inter.censor <- append(inter.censor,tem)
  volclen[jk] <- len.volc[jk]
  i <- i+1
  repi <- append(repi,rep(i,volclen[jk]))
  rep01 <- append(rep01,c(rep(0,volclen[jk]-1),1))
  Interval <- append(Interval,c(rep(max.inter,volclen[jk]-1),tem2))
  N2tem[i] <- volclen[jk]
}

dat <- cbind(inter.censor,repi,rep01,Interval)

isCensor <- dat[,3]==1
N2 <- cumsum(N2tem) # The indices for the last data points from all the volcano in dat
N1 <- c(1,N2[-length(N2)]+1) # The indices for the first data points from all the volcano in dat
Mvol <- length(Nind) # number of analogue volcanoes


# data contains the interevent times and the censored time of each volcano
# and uses a column of 1,2,3, ... to represent times from different volcanoes.

#######################################
# R2jags
library(R2jags)
library(lattice)


## Model 0: no Z. Assume all volcanic eruption records following the same renewal process
# model 0
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.45, "beta" = 1)
}

params <- c("alpha", "beta", "xi[6]")

m0 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="ModelM0Weib.jags")

m0.mcmc <- as.mcmc(m0)
summary(m0.mcmc)

m3.x14mcmc <- m0.mcmc[,4]
xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975))*1000+1500)


save(m0,file="TgriroGVPV3Weibm0-1010000i10000b.image")
save(m0.mcmc,file="TgriroGVPV3Weibm0mc-1010000i10000b.image")






