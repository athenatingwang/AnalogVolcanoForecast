## VEI 3+
source("RprogramsNatGeo/ReadGVPdata.R")

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
source("RprogramsNatGeo/denslines.R")


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
#set.seed(123)
Sys.time()
m0 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod0WeibsepZuse.jags")
Sys.time()
print(m0)
#pD = 3.9 and DIC = 2623.6
m0.mcmc <- as.mcmc(m0)
summary(m0.mcmc)

gelman.diag(m0.mcmc)

m3.x14mcmc <- m0.mcmc[,4]
xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975))*1000+1500)
  2.5%    25%    50%    75%  97.5% 
  2098   3249   6299  19100 925837 


save(m0,file="TgriroGVPV3Weibm0-1010000i10000b.image")
save(m0.mcmc,file="TgriroGVPV3Weibm0mc-1010000i10000b.image")

#get(load("Results/GVPM0-M3-EmpAnalogV3/TgriroGVPV3Weibm0-1010000i10000b.image"))
#get(load("Results/GVPM0-M3-EmpAnalogV3/TgriroGVPV3Weibm0mc-1010000i10000b.image"))


postscript("GVPV3-WeibM0sepZ-1.eps")
par(mfrow=c(2,2))
plot(m0.mcmc, trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()



## Forecasts

m3.x14mcmc <- m0.mcmc[,4]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()




## Residual analysis - Bayesian posterior

#### KS test statistic
len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam0t <- function(volci,param){ # volci: volcano number
                                    # IEk: number of eruptions from 1 to N-1
  alpha <- param[1]
  beta <- param[2]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
#  sum(-log(1-plnorm(IEtime[1:IEk],mu,sig))) # calculate Lambda_{t_k}
#  IEsort <- seq(0.001,max(IEtime),length.out=200)
#  IEsort <- c(seq(0.001,median(IEtime),length.out=len.out),seq(median(IEtime),max(IEtime),length.out=len.out))
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=alpha,scale=beta)  
  return(ksstat)
}

tem1 <- dat[volc==1,1] # Interevent times for volcano number volci
IEt <- tem1[-length(tem1)] # getting rid of the censored time
ltind <- sort(IEt);
lt <- length(ltind); 

len.chain <- length((m0.mcmc[,1])[[1]])
Lam1 <- matrix(NA,len.chain,lt1)

for (i in 1:len.chain){
  tem <- as.vector(m0.mcmc[i,][[1]])
  est1 <- tem[c(1,2)]
  Lam1[i,] <- Lam0t(1,est1)
}

lam1CI <- t(apply(Lam1,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

postscript("GVPV3TgrAlone-M0-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,cex.axis=1.4)
abline(h=0)
dev.off()





