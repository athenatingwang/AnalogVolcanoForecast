source("RprogramsNatGeo/ReadGVPdata.R")

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
  
#ind <- c(3,4,8,15)+1 # overlap M3 model greater than 0.1
# > floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975))*1000+1500)
#  0.5% 2.5%  5%   25%   50%   75%   95%  97.5%  99.5%
# 2026  2055 2095  2480  3261  4984 10930 14663  26853

ind <- c(3,4,8,15,17)+1 # overlap M3 model greater than 0.2
#> floor(quantile((xx[,1,drop=TRUE]),c(0.05,0.25,0.5,0.75,0.95))*1000+1500)
# 0.5% 2.5%  5%  25%  50%  75%  95%   97.5% 99.5%
# 2024 2048 2078 2369 2975 4242 8715  11039 19962

#ind <- c(3,4,8,9,10,14,15,17)+1 # overlap M1 and M3 model greater than 0.2
#> floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995))*1000+1500)
# 0.5%  2.5%    5%   25%   50%   75%   95% 97.5% 99.5% 
# 2023  2042  2069  2338  2901  4086  8245 10666 17976 

#ind <- c(3,4,8,9,14,15,17)+1 # overlap M1 and M3 model greater than 0.15
#> floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995))*1000+1500)
# 0.5%  2.5%    5%   25%   50%   75%   95% 97.5% 99.5% 
# 2023  2045  2071  2341  2872  4019  7706 10043 18405 


aname <- aname[-ind]

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
source("RprogramsNatGeo/denslines.R")

# model 1
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.7, "beta" = 0.3,"sigmaZ" = 1, "Z"=rep(1,Mvol))
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[6]")

s1 <- Sys.time()
m1 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod1WeibsepZuse.jags")
s2 <- Sys.time()

m1.mcmc <- as.mcmc(m1)
summary(m1.mcmc)

gelman.diag(m1.mcmc)


save(m1,file="TgriroGVPStatsV3m1-110000i10000b.image")
save(m1.mcmc,file="TgriroGVPStatsV3m1mc-110000i10000b.image")

#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m1-110000i10000b.image"))
#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m1mc-110000i10000b.image"))

postscript("GVP-StatsV3M1sepZ.eps")
par(mfrow=c(5,5))
plot(m1.mcmc[,c(1:6,16:23,7:15)], trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()

m1.Zmcmc <- m1.mcmc[,c(6,16:23,7:15)]
m1.sZmcmc <- m1.mcmc[,4]
temp <- rainbow(Mvol)
colsN <- temp 
colsN[1] <- 1
colsN[3] <- grey(0.4)

postscript("GVP-StatsV3M1postdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2), mar=c(5.1, 5.1, 1, 1))
densplot(m1.Zmcmc[,1],xlim=c(0,4),ylim=c(0,3.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
  denslines(m1.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:ceiling(Mvol/2))
  legend(1.8,3.5-i*0.3,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
for (i in (ceiling(Mvol/2)+1):Mvol)
 legend(2.8,3.5-(i-ceiling(Mvol/2))*0.3,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m1.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
dev.off()


## p value for similarities between two catalogues using models with given initials
# Chain 1
len.chain <- length((m1.Zmcmc[,1])[[1]])
LL1c1 <- matrix(NA,Mvol-1,len.chain)
tem1c1 <- p1c1 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL1c1[j-1,i] <- sum((m1.Zmcmc[i,1])[[1]]>(m1.Zmcmc[,j])[[1]])
  }
  tem1c1[j-1] <- sum(LL1c1[j-1,])/len.chain^2
  if (tem1c1[j-1]>0.5) tem1c1[j-1] <- 1-tem1c1[j-1]
  p1c1[j-1] <- 2*tem1c1[j-1]
}
# Chain 2
LL1c2 <- matrix(NA,Mvol-1,len.chain)
tem1c2 <- p1c2 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL1c2[j-1,i] <- sum((m1.Zmcmc[i,1])[[2]]>(m1.Zmcmc[,j])[[2]])
  }
  tem1c2[j-1] <- sum(LL1c2[j-1,])/len.chain^2
  if (tem1c2[j-1]>0.5) tem1c2[j-1] <- 1-tem1c2[j-1]
  p1c2[j-1] <- 2*tem1c2[j-1]
}
# Chain 3
LL1c3 <- matrix(NA,Mvol-1,len.chain)
tem1c3 <- p1c3 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL1c3[j-1,i] <- sum((m1.Zmcmc[i,1])[[3]]>(m1.Zmcmc[,j])[[3]])
  }
  tem1c3[j-1] <- sum(LL1c3[j-1,])/len.chain^2
  if (tem1c3[j-1]>0.5) tem1c3[j-1] <- 1-tem1c3[j-1]
  p1c3[j-1] <- 2*tem1c3[j-1]
}

options("scipen"=100, "digits"=4)
apply(cbind(p1c1,p1c2,p1c3),1,mean)

 [1] 0.0807822 0.0010165 0.0258891 0.1403842 0.0008079 0.0004095 0.1648379
 [8] 0.2807484 0.0297396 0.2868765 0.2329621 0.0791466 0.0300237 0.1041481
[15] 0.0234528 0.0171441 0.4327147



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

s1 <- Sys.time()
m2 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod2WeibsepZuse.jags")
s2 <- Sys.time()
s2-s1

m2.mcmc <- as.mcmc(m2)
summary(m2.mcmc)

gelman.diag(m2.mcmc)

save(m2,file="TgriroGVPStatsV3m2-110000i10000b.image")
save(m2.mcmc,file="TgriroGVPStatsV3m2mc-110000i10000b.image")

#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m2-110000i10000b.image"))
#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m2mc-110000i10000b.image"))


postscript("GVP-StatsV3M2sepZ.eps")
par(mfrow=c(5,5))
plot(m2.mcmc[,c(1:6,16:23,7:15)], trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()


m2.Zmcmc <- m2.mcmc[,c(6,16:23,7:15)]
m2.sZmcmc <- m2.mcmc[,4]
temp <- rainbow(Mvol)
colsN <- temp 
colsN[1] <- 1
colsN[3] <- grey(0.4)

postscript("GVP-StatsV3M2postdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2), mar=c(5.1, 5.1, 1, 1))
densplot(m2.Zmcmc[,1],xlim=c(0,3.5),ylim=c(0,4.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m2.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:ceiling(Mvol/2))
 legend(1.9,4.2-i*0.29,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
for (i in (ceiling(Mvol/2)+1):Mvol)
 legend(2.7,4.2-(i-ceiling(Mvol/2))*0.29,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m2.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
dev.off()



## p value for similarities between two catalogues using models with given initials
# Chain 1
len.chain <- length((m2.Zmcmc[,1])[[1]])
LL2c1 <- matrix(NA,Mvol-1,len.chain)
tem2c1 <- p2c1 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL2c1[j-1,i] <- sum((m2.Zmcmc[i,1])[[1]]>(m2.Zmcmc[,j])[[1]])
  }
  tem2c1[j-1] <- sum(LL2c1[j-1,])/len.chain^2
  if (tem2c1[j-1]>0.5) tem2c1[j-1] <- 1-tem2c1[j-1]
  p2c1[j-1] <- 2*tem2c1[j-1]
}
# Chain 2
LL2c2 <- matrix(NA,Mvol-1,len.chain)
tem2c2 <- p2c2 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL2c2[j-1,i] <- sum((m2.Zmcmc[i,1])[[2]]>(m2.Zmcmc[,j])[[2]])
  }
  tem2c2[j-1] <- sum(LL2c2[j-1,])/len.chain^2
  if (tem2c2[j-1]>0.5) tem2c2[j-1] <- 1-tem2c2[j-1]
  p2c2[j-1] <- 2*tem2c2[j-1]
}
# Chain 3
LL2c3 <- matrix(NA,Mvol-1,len.chain)
tem2c3 <- p2c3 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL2c3[j-1,i] <- sum((m2.Zmcmc[i,1])[[3]]>(m2.Zmcmc[,j])[[3]])
  }
  tem2c3[j-1] <- sum(LL2c3[j-1,])/len.chain^2
  if (tem2c3[j-1]>0.5) tem2c3[j-1] <- 1-tem2c3[j-1]
  p2c3[j-1] <- 2*tem2c3[j-1]
}

apply(cbind(p2c1,p2c2,p2c3),1,mean)

 [1] 0.8867 0.7588 0.9596 0.8187 0.2944 0.2004 0.3971 0.7579 0.9577 0.7754
[11] 0.3358 0.6286 0.7719 0.6751 0.9301 0.3624 0.4421



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

s1 <- Sys.time()
m3 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod3WeibsepZuse.jags")
s2 <- Sys.time()
s2-s1

m3.mcmc <- as.mcmc(m3)
summary(m3.mcmc)

gelman.diag(m3.mcmc)


m3.x14mcmc <- m3.mcmc[,6]
m3.x14mcmc <- m2.mcmc[,5]

xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995))*1000+1500)


library("MCMCvis")
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)


save(m3,file="TgriroGVPStatsV3m3-110000i10000b.image")
save(m3.mcmc,file="TgriroGVPStatsV3m3mc-110000i10000b.image")

#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m3-110000i10000b.image"))
#get(load("Results/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m3mc-110000i10000b.image"))


postscript("GVP-StatsV3M3sepZ-1.eps")
par(mfrow=c(5,5))
plot(m3.mcmc[,c(1:6,25,35:42,26:34)], trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()
postscript("GVP-StatsV3M3sepZ-2.eps")
par(mfrow=c(4,5))
plot(m3.mcmc[,c(7,17:24,8:16)], trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()


#m3.Ymcmc <- m3.mcmc[,c(7,18:25,8:17)]
#m3.Zmcmc <- m3.mcmc[,c(26,37:44,27:36)]

m3.Ymcmc <- m3.mcmc[,c(7,17:24,8:16)]
m3.Zmcmc <- m3.mcmc[,c(25,35:42,26:34)]
m3.sYmcmc <- m3.mcmc[,4]
m3.sZmcmc <- m3.mcmc[,5]
temp <- rainbow(Mvol)
colsN <- temp 
colsN[1] <- 1

postscript("GVP-StatsV3M3postdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5*2,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(2,2), mar=c(5.1, 5.1, 1, 1))
densplot(m3.Zmcmc[,1],xlim=c(0,3.5),ylim=c(0,6.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:ceiling(Mvol/2))
 legend(1.9,5.5-i*0.33,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
for (i in (ceiling(Mvol/2)+1):Mvol)
 legend(2.7,5.5-(i-ceiling(Mvol/2))*0.33,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m3.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)

densplot(m3.Ymcmc[,1],xlim=c(0,5.5),ylim=c(0,5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:ceiling(Mvol/2))
 legend(2.4,5-i*0.32,paste("Y",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
for (i in (ceiling(Mvol/2)+1):Mvol)
 legend(3.6,5-(i-ceiling(Mvol/2))*0.32,paste("Y",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m3.sYmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
dev.off()




## p value for similarities between two catalogues using models with given initials
# Chain 1
len.chain <- length((m3.Zmcmc[,1])[[1]])
LL3c1 <- matrix(NA,Mvol-1,len.chain)
tem3c1 <- p3c1 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL3c1[j-1,i] <- sum((m3.Zmcmc[i,1])[[1]]>(m3.Zmcmc[,j])[[1]])
  }
  tem3c1[j-1] <- sum(LL3c1[j-1,])/len.chain^2
  if (tem3c1[j-1]>0.5) tem3c1[j-1] <- 1-tem3c1[j-1]
  p3c1[j-1] <- 2*tem3c1[j-1]
}
# Chain 2
LL3c2 <- matrix(NA,Mvol-1,len.chain)
tem3c2 <- p3c2 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL3c2[j-1,i] <- sum((m3.Zmcmc[i,1])[[2]]>(m3.Zmcmc[,j])[[2]])
  }
  tem3c2[j-1] <- sum(LL3c2[j-1,])/len.chain^2
  if (tem3c2[j-1]>0.5) tem3c2[j-1] <- 1-tem3c2[j-1]
  p3c2[j-1] <- 2*tem3c2[j-1]
}
# Chain 3
LL3c3 <- matrix(NA,Mvol-1,len.chain)
tem3c3 <- p3c3 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LL3c3[j-1,i] <- sum((m3.Zmcmc[i,1])[[3]]>(m3.Zmcmc[,j])[[3]])
  }
  tem3c3[j-1] <- sum(LL3c3[j-1,])/len.chain^2
  if (tem3c3[j-1]>0.5) tem3c3[j-1] <- 1-tem3c3[j-1]
  p3c3[j-1] <- 2*tem3c3[j-1]
}

apply(cbind(p3c1,p3c2,p3c3),1,mean)


 

## p value for similarities between two catalogues (Y) using models with given initials
# Chain 1
LLY3c1 <- matrix(NA,Mvol-1,len.chain)
temY3c1 <- pY3c1 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LLY3c1[j-1,i] <- sum((m3.Ymcmc[i,1])[[1]]>(m3.Ymcmc[,j])[[1]])
  }
  temY3c1[j-1] <- sum(LLY3c1[j-1,])/len.chain^2
  if (temY3c1[j-1]>0.5) temY3c1[j-1] <- 1-temY3c1[j-1]
  pY3c1[j-1] <- 2*temY3c1[j-1]
}
# Chain 2
LLY3c2 <- matrix(NA,Mvol-1,len.chain)
temY3c2 <- pY3c2 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LLY3c2[j-1,i] <- sum((m3.Ymcmc[i,1])[[2]]>(m3.Ymcmc[,j])[[2]])
  }
  temY3c2[j-1] <- sum(LLY3c2[j-1,])/len.chain^2
  if (temY3c2[j-1]>0.5) temY3c2[j-1] <- 1-temY3c2[j-1]
  pY3c2[j-1] <- 2*temY3c2[j-1]
}
# Chain 3
LLY3c3 <- matrix(NA,Mvol-1,len.chain)
temY3c3 <- pY3c3 <- NULL
for (j in 2:Mvol){
  for (i in 1:len.chain){
    LLY3c3[j-1,i] <- sum((m3.Ymcmc[i,1])[[3]]>(m3.Ymcmc[,j])[[3]])
  }
  temY3c3[j-1] <- sum(LLY3c3[j-1,])/len.chain^2
  if (temY3c3[j-1]>0.5) temY3c3[j-1] <- 1-temY3c3[j-1]
  pY3c3[j-1] <- 2*temY3c3[j-1]
}

apply(cbind(pY3c1,pY3c2,pY3c3),1,mean)




## Residual analysis - Bayesian posterior

#### KS test statistic
#len.out <- 300
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





gelman.diag(m1.mcmc)
gelman.diag(m2.mcmc)
gelman.diag(m3.mcmc)


options("scipen"=100, "digits"=4)
apply(cbind(p1c1,p1c2,p1c3),1,mean)
apply(cbind(p2c1,p2c2,p2c3),1,mean)
apply(cbind(p3c1,p3c2,p3c3),1,mean)
apply(cbind(pY3c1,pY3c2,pY3c3),1,mean)






## Forecast x[6]
postscript("GVPV3-Statsrec-M123-Forecast.eps",paper="special",
    width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,4.5,2.2,1))
m1.x14mcmc <- m1.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995)),
     lty="dashed",lwd=2,col="blue")
#quantile((xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995))*1000+1798
# 2022.378  2347.619  3121.109  5432.091 58934.643 
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.5,0.995)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.5,0.995))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(3.8,0.42,"(a)",cex=1.5)

m2.x14mcmc <- m2.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995)),
     lty="dashed",lwd=2,col="blue")
#quantile((xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995))*1000+1798
#  2024.468   2585.238   4247.916  10331.148 231214.742 
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.5,0.995)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.5,0.995))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(2.2,1.1,"(b)",cex=1.5)

m3.x14mcmc <- m3.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995)),
     lty="dashed",lwd=2,col="blue")
#quantile((xx[,1,drop=TRUE]),c(0.005,0.25,0.5,0.75,0.995))*1000+1798
# 2022.287   2373.197   3257.082   6203.279 136560.033
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.005,0.5,0.995)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.5,0.995))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(5.8,0.37,"(c)",cex=1.5)
dev.off()






