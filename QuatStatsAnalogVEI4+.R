## VEI 4+
tgr3 <- read.csv("QuaternaryRecords/Tongariro.csv",header=TRUE)
tgr3$VEI[is.na(tgr3$VEI)] <- 4.5
tgr <- tgr3[tgr3$VEI>=3.99,]
ths3 <- read.csv("QuaternaryRecords/ThreeSisters.csv",header=TRUE)
ths3$VEI[is.na(ths3$VEI)] <- 4.5
ths <- ths3[ths3$VEI>=3.99,]
pey3  <- read.csv("QuaternaryRecords/Puyehue-CordonCaulle.csv",header=TRUE)
pey3$VEI[is.na(pey3$VEI)] <- 4.5
pey <- pey3[pey3$VEI>=3.99,]


tgrIE <- -diff(tgr$AgeyrBP)/1000
thsIE <- -diff(ths$AgeyrBP)/1000
peyIE <- -diff(pey$AgeyrBP)/1000

ytgr <- -(1950-2018-tgr$AgeyrBP[length(tgr$AgeyrBP)])/1000
yths <- -(1950-2018-ths$AgeyrBP[length(ths$AgeyrBP)])/1000
ypey <- -(1950-2018-pey$AgeyrBP[length(pey$AgeyrBP)])/1000

temp <- c(tgrIE,ytgr,thsIE,yths,peyIE,ypey)
tmtgr <- max(tgrIE)+1
tmths <- max(thsIE)+1
tmpey <- max(peyIE)+1
dat <- cbind(c(tgrIE,NA,thsIE,NA,peyIE,NA),
       c(rep(1,length(tgr$AgeyrBP)),rep(2,length(ths$AgeyrBP)),
          rep(3,length(pey$AgeyrBP))),
          c(rep(0,length(tgr$AgeyrBP)-1),1,rep(0,length(ths$AgeyrBP)-1),1,
          rep(0,length(pey$AgeyrBP)-1),1),
		  c(rep(tmtgr,length(tgr$AgeyrBP)-1),ytgr,rep(tmths,length(ths$AgeyrBP)-1),
		  yths,rep(tmpey,length(pey$AgeyrBP)-1),ypey))
isCensor <- dat[,3]==1
N1 <- N2 <- NULL
N2[1] <- length(tgr$AgeyrBP)
N2[2] <- length(tgr$AgeyrBP)+length(ths$AgeyrBP)
N2[3] <- length(dat[,1])
N1[1] <- 1
N1[2] <- N2[1]+1
N1[3] <- N2[2]+1
Mvol <- 3

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
  list("alpha" = 0.45, "beta" = 2.5, "sigmaZ" = 1, "Z"=rep(1,Mvol))#,"xi"=xiInit)
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[9]")

m1 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="ModelM1Weib.jags")

m1.mcmc <- as.mcmc(m1)
summary(m1.mcmc)

save(m1,file="TgriroLongRecWeibm1-1010000i10000b.image")
save(m1.mcmc,file="TgriroLongRecWeibm1mc-1010000i10000b.image")

m1.Zmcmc <- m1.mcmc[,c(6:8)]
m1.sZmcmc <- m1.mcmc[,4]


## Residual analysis - Bayesian posterior

#### KS test statistic
len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam1t <- function(volci,param){ # volci: volcano number
                                    # IEk: number of eruptions from 1 to N-1
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=alpha,scale=Z^(-1/alpha)*beta)  
  return(ksstat)
}

ltind1 <- sort(tgrIE);
ltind2 <- sort(thsIE);
ltind3 <- sort(peyIE)
lt1 <- length(ltind1); 
lt2 <- length(ltind2); 
lt3 <- length(ltind3); 

len.chain <- length((m1.Zmcmc[,1])[[1]])
Lam1 <- matrix(NA,len.chain,lt1)
Lam2 <- matrix(NA,len.chain,lt2); Lam3 <- matrix(NA,len.chain,lt3) 

for (i in 1:len.chain){
  tem <- as.vector(m1.mcmc[i,][[1]])
  est1 <- tem[c(1,2,6)]
  est2 <- tem[c(1,2,7)]
  est3 <- tem[c(1,2,8)]
  Lam1[i,] <- Lam1t(1,est1)
  Lam2[i,] <- Lam1t(2,est2)
  Lam3[i,] <- Lam1t(3,est3)
}

lam1CI <- t(apply(Lam1,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam2CI <- t(apply(Lam2,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam3CI <- t(apply(Lam3,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

postscript("LongRecAnalogWeib-M1-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),
        c(0.001,0.01,0.1,1,5),cex.axis=1.4)
text(-4.8,0.59,"(a)",cex=1.5)

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.01,0.1,1,10,100)),
        c(0.01,0.1,1,10,100),cex.axis=1.4)
text(-2,0.23,"(b)",cex=1.5)

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),
        c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
text(-7,0.19,"(c)",cex=1.5)
dev.off()




# model 2
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.45, "beta" = 2.5, "sigmaZ" = 0.5, "Z"=rep(1,Mvol))#, "xi"=xiInit)
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[9]")

m2 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="ModelM2Weib.jags")

m2.mcmc <- as.mcmc(m2)
summary(m2.mcmc)

save(m2,file="TgriroLongRecWeibm2-1010000i10000b.image")
save(m2.mcmc,file="TgriroLongRecWeibm2mc-1010000i10000b.image")

m2.Zmcmc <- m2.mcmc[,c(6:8)]
m2.sZmcmc <- m2.mcmc[,4]


## Residual analysis - Bayesian posterior

#### KS test statistic
len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam1t <- function(volci,param){ # volci: volcano number
                                    # IEk: number of eruptions from 1 to N-1
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  temp <- dat[volc==volci,1] # Interevent times for volcano number volci
  IEtime <- temp[-length(temp)] # getting rid of the censored time
  IEsort <- sort(IEtime)
  ksstat <- ecdf(IEtime)(IEsort)-pweibull(IEsort,shape=Z*alpha,scale=Z^(1/Z/alpha)*beta^(1/Z))  
  return(ksstat)
}

ltind1 <- sort(tgrIE);
ltind2 <- sort(thsIE);
ltind3 <- sort(peyIE)
lt1 <- length(ltind1); 
lt2 <- length(ltind2); 
lt3 <- length(ltind3); 

len.chain <- length((m2.Zmcmc[,1])[[1]])
Lam1 <- matrix(NA,len.chain,lt1)
Lam2 <- matrix(NA,len.chain,lt2); Lam3 <- matrix(NA,len.chain,lt3) 

for (i in 1:len.chain){
  tem <- as.vector(m2.mcmc[i,][[1]])
  est1 <- tem[c(1,2,6)]
  est2 <- tem[c(1,2,7)]
  est3 <- tem[c(1,2,8)]
  Lam1[i,] <- Lam1t(1,est1)
  Lam2[i,] <- Lam1t(2,est2)
  Lam3[i,] <- Lam1t(3,est3)
}

lam1CI <- t(apply(Lam1,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam2CI <- t(apply(Lam2,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam3CI <- t(apply(Lam3,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

postscript("LongRecAnalogWeib-M2-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),
        c(0.001,0.01,0.1,1,5),cex.axis=1.4)
text(-4.8,0.59,"(d)",cex=1.5)

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.01,0.1,1,10,100)),
        c(0.01,0.1,1,10,100),cex.axis=1.4)
text(-2,0.21,"(e)",cex=1.5)

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),
        c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
text(-7,0.2,"(f)",cex=1.5)
dev.off()




# model 3
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.45, "beta" = 2.5, "sigmaZ" = 1, "Z"=rep(1,Mvol),
       "sigmaY" = 1, "Y"=rep(1,Mvol))#, 
#       "xi"=xiInit)
}

params <- c("alpha", "beta", "sigmaZ", "Z", "sigmaY", "Y","xi[9]")

m3 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="ModelM3Weib.jags")

m3.mcmc <- as.mcmc(m3)
summary(m3.mcmc)

save(m3,file="TgriroLongRecWeibm3-1010000i10000b.image")
save(m3.mcmc,file="TgriroLongRecWeibm3mc-1010000i10000b.image")


m3.Ymcmc <- m3.mcmc[,c(7:9)]
m3.Zmcmc <- m3.mcmc[,c(10:12)]
m3.sYmcmc <- m3.mcmc[,4]
m3.sZmcmc <- m3.mcmc[,5]


## Residual analysis - Bayesian posterior

#### KS test statistic
len.out <- 300
volc <- dat[,2] # A vector indicating which volcanoes 
## (as a number ordered according to volcano name alphabetacally)
Lam1t <- function(volci,param){ # volci: volcano number
                                    # IEk: number of eruptions from 1 to N-1
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

ltind1 <- sort(tgrIE);
ltind2 <- sort(thsIE);
ltind3 <- sort(peyIE)
lt1 <- length(ltind1); 
lt2 <- length(ltind2); 
lt3 <- length(ltind3); 

len.chain <- length((m3.Zmcmc[,1])[[1]])
Lam1 <- matrix(NA,len.chain,lt1)
Lam2 <- matrix(NA,len.chain,lt2); Lam3 <- matrix(NA,len.chain,lt3) 

for (i in 1:len.chain){
  tem <- as.vector(m3.mcmc[i,][[1]])
  est1 <- tem[c(1,2,10,7)]
  est2 <- tem[c(1,2,11,8)]
  est3 <- tem[c(1,2,12,9)]
  Lam1[i,] <- Lam1t(1,est1)
  Lam2[i,] <- Lam1t(2,est2)
  Lam3[i,] <- Lam1t(3,est3)
}

lam1CI <- t(apply(Lam1,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam2CI <- t(apply(Lam2,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))
lam3CI <- t(apply(Lam3,2,quantile,c(0.005,0.25,0.5,0.75,0.995)))

postscript("LongRecAnalogWeib-M3-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),
        c(0.001,0.01,0.1,1,5),cex.axis=1.4)
text(-4.8,0.59,"(g)",cex=1.5)

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.01,0.1,1,10,100)),
        c(0.01,0.1,1,10,100),cex.axis=1.4)
text(-2,0.255,"(h)",cex=1.5)

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,xaxt="n")
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),
        c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
text(-7,0.19,"(i)",cex=1.5)
dev.off()


