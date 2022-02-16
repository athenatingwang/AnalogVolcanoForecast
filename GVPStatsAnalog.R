source("ReadGVPdata.R")

## Using eruption records with minimum VEI 3

# From analysis of all empirical VEI 3+ records, Number c(3,4,8,15)+1 volcanoes 
# are not analogue volcanoes of Tongariro. We therefore remove them from the 
# analogue set and redo the analysis.

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

ind <- c(3,4,8,15,17)+1 

aname <- aname[-ind]  # This subset of volcanoes satisfies that the proportion that the MCMC samples of the posterior of Z_i overlap with 
                      # Tongariro values (i.e., Z_1) for model M3 is greater than 0.2.

Nind1 <- c(21,1:20,22:23)
Nind <- Nind1[-ind]

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
N2tem <- 0
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

# model 1
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.7, "beta" = 0.3,"sigmaZ" = 1, "Z"=rep(1,Mvol))
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[6]")

m1 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="ModelM1Weib.jags")

m1.mcmc <- as.mcmc(m1)
summary(m1.mcmc)

save(m1,file="TgriroGVPStatsV3m1-110000i10000b.image")
save(m1.mcmc,file="TgriroGVPStatsV3m1mc-110000i10000b.image")

m1.Zmcmc <- m1.mcmc[,c(6,16:23,7:15)]
m1.sZmcmc <- m1.mcmc[,4]

## Residual analysis - Bayesian posterior

#### KS test statistic
#len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam1t <- function(volci,param){ # volci: volcano number
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=alpha,scale=Z^(-1/alpha)*beta)  
  return(ksstat)
}

len.chain <- length((m1.Zmcmc[,1])[[1]])
postscript("GVP-StatsV3-M1-Resid-log.eps",paper="special",width=2*5*cos (35.4/180*pi)/0.612,
    height=0.5*5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(5,4),mar=c(4.5,5,1,1))

for (jk in 1:length(Nind)){
  tem1 <- dat[volc==jk,1] # Interevent times for volcano number volci
  IEt <- tem1[-length(tem1)] # getting rid of the censored time
  ltind <- sort(IEt);
  lt <- length(ltind); 
  Lam <- matrix(NA,len.chain,lt)
  for (i in 1:len.chain){
    temab <- as.vector(m1.mcmc[i,][[1]])
    temZ <- as.vector(m1.Zmcmc[i,][[1]])
    est1 <- c(temab[c(1,2)],temZ[jk])
    Lam[i,] <- Lam1t(jk,est1)
  }
  lamCI <- t(apply(Lam,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

  matplot(log(ltind[-lt]),lamCI[-lt,],ylim=range(lamCI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    main=paste(aname[jk],"(GVP)",sep=" "),
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100)),
        c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100),cex.axis=1.4)
}
dev.off()




# model 2		
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.7, "beta" = 0.3,"sigmaZ" = 1, "Z"=rep(1,Mvol))
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[6]")

m2 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="ModelM2Weib.jags")

m2.mcmc <- as.mcmc(m2)
summary(m2.mcmc)

save(m2,file="TgriroGVPStatsV3m2-110000i10000b.image")
save(m2.mcmc,file="TgriroGVPStatsV3m2mc-110000i10000b.image")

m2.Zmcmc <- m2.mcmc[,c(6,16:23,7:15)]
m2.sZmcmc <- m2.mcmc[,4]

## Residual analysis - Bayesian posterior

#### KS test statistic
#len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam2t <- function(volci,param){ # volci: volcano number
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=Z*alpha,scale=Z^(1/Z/alpha)*beta^(1/Z))  
  return(ksstat)
}

len.chain <- length((m2.Zmcmc[,1])[[1]])
postscript("GVP-StatsV3-M2-Resid-log.eps",paper="special",width=2*5*cos (35.4/180*pi)/0.612,
    height=0.5*5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(5,4),mar=c(4.5,5,1,1))

for (jk in 1:length(Nind)){
  tem1 <- dat[volc==jk,1] # Interevent times for volcano number volci
  IEt <- tem1[-length(tem1)] # getting rid of the censored time
  ltind <- sort(IEt);
  lt <- length(ltind); 
  Lam <- matrix(NA,len.chain,lt)
  for (i in 1:len.chain){
    temab <- as.vector(m2.mcmc[i,][[1]])
    temZ <- as.vector(m2.Zmcmc[i,][[1]])
    est1 <- c(temab[c(1,2)],temZ[jk])
    Lam[i,] <- Lam2t(jk,est1)
  }
  lamCI <- t(apply(Lam,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

  matplot(log(ltind[-lt]),lamCI[-lt,],ylim=range(lamCI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    main=paste(aname[jk],"(GVP)",sep=" "),
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
 abline(h=0)
 axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100)),
         c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100),cex.axis=1.4)
}
dev.off()




# model 3
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.8, "beta" = 0.3, "sigmaZ" = 1, "Z"=rep(1,Mvol), 
       "sigmaY" = 1, "Y"=rep(1,Mvol))
}

params <- c("alpha", "beta", "sigmaY", "Y", "sigmaZ", "Z","xi[6]")

m3 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="ModelM3Weib.jags")

m3.mcmc <- as.mcmc(m3)
summary(m3.mcmc)

library("MCMCvis")
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)

save(m3,file="TgriroGVPStatsV3m3-110000i10000b.image")
save(m3.mcmc,file="TgriroGVPStatsV3m3mc-110000i10000b.image")


m3.Ymcmc <- m3.mcmc[,c(7,17:24,8:16)]
m3.Zmcmc <- m3.mcmc[,c(25,35:42,26:34)]
m3.sYmcmc <- m3.mcmc[,4]
m3.sZmcmc <- m3.mcmc[,5]


## Residual analysis - Bayesian posterior

#### KS test statistic
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam3t <- function(volci,param){ # volci: volcano number
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  Y <- param[4]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=Z*alpha,scale=Z^(1/Z/alpha)*Y^(-1/Z/alpha)*beta^(1/Z))  
  return(ksstat)
}

len.chain <- length((m3.Zmcmc[,1])[[1]])
postscript("GVP-StatsV3-M3-Resid-log.eps",paper="special",width=2*5*cos (35.4/180*pi)/0.612,
    height=0.5*5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(5,4),mar=c(4.5,5,1,1))

a2h <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)")
kk <- 1
for (jk in 1:length(Nind)){
  tem1 <- dat[volc==jk,1] # Interevent times for volcano number volci
  IEt <- tem1[-length(tem1)] # getting rid of the censored time
  ltind <- sort(IEt);
  lt <- length(ltind); 
  Lam <- matrix(NA,len.chain,lt)
  for (i in 1:len.chain){
    temab <- as.vector(m3.mcmc[i,][[1]])
    temZ <- as.vector(m3.Zmcmc[i,][[1]])
    temY <- as.vector(m3.Ymcmc[i,][[1]])
    est1 <- c(temab[c(1,2)],temZ[jk],temY[jk])
    Lam[i,] <- Lam3t(jk,est1)
  }
  lamCI <- t(apply(Lam,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

  matplot(log(ltind[-lt]),lamCI[-lt,],ylim=range(lamCI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    main=paste(aname[jk],"(GVP)",sep=" "),
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
  abline(h=0)
  axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100)),
          c(0.00001,0.0001,0.001,0.01,0.1,0.5,1,10,100),cex.axis=1.4)
#  text(max(log(ltind))-(max(log(ltind))-min(log(ltind)))/5,
#       max(lamCI)-(max(lamCI)-min(lamCI))/10,a2h[kk],cex=1.5)
kk <- kk+1
}
dev.off()


