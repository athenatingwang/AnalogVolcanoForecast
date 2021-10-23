## VEI 3+
tgr3 <- read.csv("AnaloguesFromMark/Tongariro.csv",header=TRUE)
tgr3$VEI[is.na(tgr3$VEI)] <- 4.5
tgr <- tgr3[tgr3$VEI>=2.99,]

tgrIE <- -diff(tgr$AgeyrBP)/1000

ytgr <- -(1950-2018-tgr$AgeyrBP[length(tgr$AgeyrBP)])/1000

temp <- c(tgrIE,ytgr)
tmtgr <- max(tgrIE)+1
dat <- cbind(c(tgrIE,NA),
       rep(1,length(tgr$AgeyrBP)),
          c(rep(0,length(tgr$AgeyrBP)-1),1),
		  c(rep(tmtgr,length(tgr$AgeyrBP)-1),ytgr))
isCensor <- dat[,3]==1
N1 <- N2 <- NULL
N2 <- length(tgr$AgeyrBP)
N1 <- 1
Mvol <- 1

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

params <- c("alpha", "beta", "xi[14]")
#set.seed(123)
s1 <- Sys.time()
m0 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod0WeibsepZuse.jags")
s2 <- Sys.time()
s2-s1

m0.mcmc <- as.mcmc(m0)
summary(m0.mcmc)

gelman.diag(m0.mcmc)

m3.x14mcmc <- m0.mcmc[,4]
xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975))*1000+1798)
 2.5%   25%   50%   75% 97.5% 
 2041  2387  3329  6733 70027 



save(m0,file="TgriroLongRecV3Weibm0-1010000i10000b.image")
save(m0.mcmc,file="TgriroLongRecV3Weibm0mc-1010000i10000b.image")

#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm0-1010000i10000b.image"))
#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm0mc-1010000i10000b.image"))


postscript("LongRecV3-WeibM0sepZ-1.eps")
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
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.05,0.25,0.5,0.75,0.95)),
     lty="dashed",lwd=2,col="blue")
#quantile((xx[,1,drop=TRUE]),c(0.05,0.25,0.5,0.75,0.95))*1000+1798
#  2066.034  2373.197  3257.082  6203.279 25217.844 
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.05,0.5,0.95)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.05,0.5,0.95))*1000+1798),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1798),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1798),
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

ltind1 <- sort(tgrIE);
lt1 <- length(ltind1); 

len.chain <- length((m0.mcmc[,1])[[1]])
Lam1 <- matrix(NA,len.chain,lt1)

for (i in 1:len.chain){
  tem <- as.vector(m0.mcmc[i,][[1]])
  est1 <- tem[c(1,2)]
  Lam1[i,] <- Lam0t(1,est1)
}

lam1CI <- t(apply(Lam1,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

postscript("LongRecV3AnalogWeib-M0-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,cex.axis=1.4)
abline(h=0)
dev.off()







