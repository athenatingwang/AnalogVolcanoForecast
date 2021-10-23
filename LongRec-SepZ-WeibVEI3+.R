## Highest posterior density (HPD) credible interval
# 95% and 50%
library("MCMCvis")
# Model 1
MCMCsummary(m1.mcmc,HPD=TRUE,hpd_prob=0.25,round=3)
            mean     sd 50%_HPDL 25%_HPDU Rhat  n.eff
xi[14]     3.959  9.399    0.220    0.550 1.00 150000
MCMCsummary(m1.mcmc,HPD=TRUE,hpd_prob=0.5,round=3)
            mean     sd 50%_HPDL 50%_HPDU Rhat  n.eff
xi[14]     3.959  9.399    0.220    1.323 1.00 150000
MCMCsummary(m1.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)
            mean     sd 95%_HPDL 95%_HPDU Rhat  n.eff
xi[14]     3.959  9.399    0.220   15.540 1.00 150000

forex.m1 <- 1000*scan()
0.220    0.550  0.220    1.323 0.220   15.540

forex.m1+1798
[1]  2018  2348  2018  3121  2018 17338


# Model 2
MCMCsummary(m2.mcmc,HPD=TRUE,hpd_prob=0.25,round=3)
            mean     sd 25%_HPDL 25%_HPDU Rhat  n.eff
xi[14]    11.931 46.802    0.220    0.787 1.00 142194
MCMCsummary(m2.mcmc,HPD=TRUE,hpd_prob=0.5,round=3)
            mean     sd 50%_HPDL 50%_HPDU Rhat  n.eff
xi[14]    11.931 46.802    0.220    2.450 1.00 142194
MCMCsummary(m2.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)
            mean     sd 95%_HPDL 95%_HPDU Rhat  n.eff
xi[14]    11.931 46.802    0.220   46.775 1.00 142194

forex.m2 <- 1000*scan()
0.220    0.787  0.220    2.450  0.220   46.775

forex.m2+1798
[1]  2018  2585  2018  4248  2018 48573


# Model 3
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.25,round=3)
            mean     sd 25%_HPDL 25%_HPDU Rhat  n.eff
xi[14]     7.543 225.587    0.220    0.575 1.13 150000
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.5,round=3)
            mean     sd 50%_HPDL 50%_HPDU Rhat  n.eff
xi[14]     7.543 225.587    0.220    1.459 1.13 150000
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)
            mean     sd 95%_HPDL 95%_HPDU Rhat  n.eff
xi[14]     7.543 225.587    0.220   23.423 1.13 150000

forex.m3 <- 1000*scan()
0.220    0.575  0.220    1.459  0.220   23.423

forex.m3+1798
[1]  2018  2373  2018  3257  2018 25221




## VEI 3+
tgr3 <- read.csv("AnaloguesFromMark/Tongariro.csv",header=TRUE)
tgr3$VEI[is.na(tgr3$VEI)] <- 4.5
tgr <- tgr3[tgr3$VEI>=2.99,]
ths3 <- read.csv("AnaloguesFromMark/ThreeSisters.csv",header=TRUE)
ths3$VEI[is.na(ths3$VEI)] <- 4.5
ths <- ths3[ths3$VEI>=2.99,]
#ths <- ths[-(1:13),]
pey3  <- read.csv("AnaloguesFromMark/Puyehue-CordonCaulle.csv",header=TRUE)
pey3$VEI[is.na(pey3$VEI)] <- 4.5
pey <- pey3[pey3$VEI>=2.99,]
#pey <- pey[-(1:20),]


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

# N1: The indices for the first data points from all the volcano in dat
# N2: The indices for the last data points from all the volcano in dat
N1 <- N2 <- NULL 
N2[1] <- length(tgr$AgeyrBP)
N2[2] <- length(tgr$AgeyrBP)+length(ths$AgeyrBP)
N2[3] <- length(dat[,1])
N1[1] <- 1
N1[2] <- N2[1]+1
N1[3] <- N2[2]+1
Mvol <- 3 # number of analogue volcanoes

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
  list("alpha" = 0.45, "beta" = 2.5, "sigmaZ" = 1, "Z"=rep(1,Mvol))#,"xi"=xiInit)
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[14]")

s1 <- Sys.time()
m1 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod1WeibsepZuse.jags")
s2 <- Sys.time()

m1.mcmc <- as.mcmc(m1)
summary(m1.mcmc)

gelman.diag(m1.mcmc)

save(m1,file="TgriroLongRecV3Weibm1-1010000i10000b.image")
save(m1.mcmc,file="TgriroLongRecV3Weibm1mc-1010000i10000b.image")

#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm1-1010000i10000b.image"))
#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm1mc-1010000i10000b.image"))


postscript("LongRecV3-WeibM1sepZ-1.eps")
par(mfrow=c(2,4))
plot(m1.mcmc, trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()

m1.Zmcmc <- m1.mcmc[,c(6:8)]
m1.sZmcmc <- m1.mcmc[,4]
temp <- rainbow(Mvol)
colsN <- temp #[c(16,9,5,11,17,7,12,3 ,2 ,1,14 ,6,13,15 ,10,8,4)]
colsN[1] <- 1

postscript("LongRecV3-WeibM1sepZpostdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2), mar=c(5.1, 5.1, 1, 1))
densplot(m1.Zmcmc[,1],ylim=c(0,3),xlim=c(0,4),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:3)
denslines(m1.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:3)
 legend(3,2.8-i*0.29,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m1.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
dev.off()



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
#  sum(-log(1-plnorm(IEtime[1:IEk],mu,sig))) # calculate Lambda_{t_k}
#  IEsort <- seq(0.001,max(IEtime),length.out=200)
#  IEsort <- c(seq(0.001,median(IEtime),length.out=len.out),seq(median(IEtime),max(IEtime),length.out=len.out))
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

postscript("LongRecV3AnalogWeib-M1-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),c(0.001,0.01,0.1,1,5),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-6.6,0.5,"(a)",cex=1.5)
box()

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.1,1,10,100,150)),c(0.1,1,10,100,150),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-2,0.23,"(b)",cex=1.5)
box()

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-7,0.19,"(c)",cex=1.5)
box()
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





# model 2
data1 <- list("xi" = dat[,1], "isCensored" = dat[,3], "censLimVec"=dat[,4], 
              "N1" = N1, "N2" = N2,"Mvol"=Mvol)
xiInit <- rep( NA , length(dat[,1]) )
xiInit[isCensor] <- (dat[,4])[isCensor]+1
vinits <- function(){
  list("alpha" = 0.45, "beta" = 2.5, "sigmaZ" = 0.5, "Z"=rep(1,Mvol))#, "xi"=xiInit)
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[14]")
#set.seed(123)
Sys.time()
m2 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod2WeibsepZuse.jags")
Sys.time()
print(m2)
#plot(m2)
#traceplot(m2)
m2.mcmc <- as.mcmc(m2)
summary(m2.mcmc)

gelman.diag(m2.mcmc)


save(m2,file="TgriroLongRecV3Weibm2-1010000i10000b.image")
save(m2.mcmc,file="TgriroLongRecV3Weibm2mc-1010000i10000b.image")

#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm2-1010000i10000b.image"))
#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm2mc-1010000i10000b.image"))


postscript("LongRecV3-WeibM2sepZ-1.eps")
par(mfrow=c(2,4))
plot(m2.mcmc, trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()

m2.Zmcmc <- m2.mcmc[,c(6:8)]
m2.sZmcmc <- m2.mcmc[,4]
temp <- rainbow(Mvol)
colsN <- temp #[c(16,9,5,11,17,7,12,3 ,2 ,1,14 ,6,13,15 ,10,8,4)]
colsN[1] <- 1

postscript("LongRecV3-WeibM2sepZpostdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2), mar=c(5.1, 5.1, 1, 1))
densplot(m2.Zmcmc[,1],xlim=c(0,2),ylim=c(0,5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:3)
denslines(m2.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:3)
 legend(1.4,5-i*0.5,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m2.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
dev.off()



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
#  sum(-log(1-plnorm(IEtime[1:IEk],mu,sig))) # calculate Lambda_{t_k}
#  IEsort <- seq(0.001,max(IEtime),length.out=200)
#  IEsort <- c(seq(0.001,median(IEtime),length.out=len.out),seq(median(IEtime),max(IEtime),length.out=len.out))
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

postscript("LongRecV3AnalogWeib-M2-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),c(0.001,0.01,0.1,1,5),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-6.6,0.5,"(d)",cex=1.5)
box()

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.1,1,10,100,150)),c(0.1,1,10,100,150),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-2,0.2,"(e)",cex=1.5)
box()

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-7,0.19,"(f)",cex=1.5)
box()
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

params <- c("alpha", "beta", "sigmaZ", "Z", "sigmaY", "Y","xi[14]")
#set.seed(123)
Sys.time()
m3 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 1010000,n.burnin=10000,n.thin=20,
          model.file="RprogramsNatGeo/GVPmod3WeibsepZuse.jags")
Sys.time()
print(m3)
m3.mcmc <- as.mcmc(m3)
summary(m3.mcmc)

gelman.diag(m3.mcmc)

save(m3,file="TgriroLongRecV3Weibm3-1010000i10000b.image")
save(m3.mcmc,file="TgriroLongRecV3Weibm3mc-1010000i10000b.image")

#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm3-1010000i10000b.image"))
#get(load("Results/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm3mc-1010000i10000b.image"))


postscript("LongRecV3-WeibM3sepZ-1.eps")
par(mfrow=c(3,4))
plot(m3.mcmc, trace = TRUE, density = FALSE,auto.layout =FALSE)
dev.off()

m3.Ymcmc <- m3.mcmc[,c(7:9)]
m3.Zmcmc <- m3.mcmc[,c(10:12)]
m3.sYmcmc <- m3.mcmc[,4]
m3.sZmcmc <- m3.mcmc[,5]
temp <- rainbow(Mvol)
colsN <- temp #[c(16,9,5,11,17,7,12,3 ,2 ,1,14 ,6,13,15 ,10,8,4)]
colsN[1] <- 1

postscript("LongRecV3-WeibM3sepZpostdens.eps",paper="special",
      width=5*cos (35.4/180*pi)/0.612*2,height=5*2,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(2,2), mar=c(5.1, 5.1, 1, 1))
densplot(m3.Zmcmc[,1],xlim=c(0,1.8),ylim=c(0,5.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:3)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:3)
 legend(1.2,5-i*0.5,paste("Z",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m3.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)

densplot(m3.Ymcmc[,1],xlim=c(0,5),ylim=c(0,2.6),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:3)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
for (i in 1:3)
 legend(3.2,2.3-i*0.22,paste("Y",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.2,lwd=1.8)
densplot(m3.sYmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
dev.off()



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
#  sum(-log(1-plnorm(IEtime[1:IEk],mu,sig))) # calculate Lambda_{t_k}
#  IEsort <- seq(0.001,max(IEtime),length.out=200)
#  IEsort <- c(seq(0.001,median(IEtime),length.out=len.out),seq(median(IEtime),max(IEtime),length.out=len.out))
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

postscript("LongRecV3AnalogWeib-M3-Resid-sqrt.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,5,1,1))
matplot(log(ltind1[-lt1]),lam1CI[-lt1,],ylim=range(lam1CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n1](x)-F[1](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.001,0.01,0.1,1,5)),c(0.001,0.01,0.1,1,5),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-6.6,0.49,"(g)",cex=1.5)
box()

matplot(log(ltind2[-lt2]),lam2CI[-lt2,],ylim=range(lam2CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n2](x)-F[2](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.1,1,10,100,150)),c(0.1,1,10,100,150),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-2,0.26,"(h)",cex=1.5)
box()

matplot(log(ltind3[-lt3]),lam3CI[-lt3,],ylim=range(lam3CI),type="l",col=c(1,2,4,2,1),
    lty=c(2,3,1,3,2),xlab="Interevent times (x1000 years)",
    ylab=expression(F[n3](x)-F[3](x)),lwd=1.5,cex.lab=1.4,axes=F)
abline(h=0)
axis(1,at=log(c(0.0001,0.001,0.01,0.1,1,10,50)),c(0.0001,0.001,0.01,0.1,1,10,50),cex.axis=1.4)
axis(2,cex.axis=1.4)
text(-7,0.18,"(i)",cex=1.5)
box()
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








## Forecast x[14]
postscript("LongRecV3AnalogWeib-M123-Forecast.eps",paper="special",
    width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,4.5,2.2,1))
m1.x14mcmc <- m1.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     c(2040,3121,25889),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     c(3121),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     c(2347,5432),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(6.8,0.28,"(a)",cex=1.5)

m2.x14mcmc <- m2.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     c(2050,4247,81541),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     c(2585,10331),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(6.8,0.215,"(b)",cex=1.5)

m3.x14mcmc <- m3.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     c(2041,3257,42612),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     c(3257),
     line=-0,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     c(2373,6203),
     line=-0.2,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
box()
text(6.8,0.265,"(c)",cex=1.5)
dev.off()




gelman.diag(m1.mcmc)
gelman.diag(m2.mcmc)
gelman.diag(m3.mcmc)

options("scipen"=100, "digits"=4)
apply(cbind(p1c1,p1c2,p1c3),1,mean)
apply(cbind(p2c1,p2c2,p2c3),1,mean)
apply(cbind(p3c1,p3c2,p3c3),1,mean)
apply(cbind(pY3c1,pY3c2,pY3c3),1,mean)













## Bayesian predictive distributions for the waiting time of the next eruption 
## xifore.m1 is the waiting time of the forecasted event from model 1
## xifore.m2 is the waiting time of the forecasted event from model 2
## xifore.m3 is the waiting time of the forecasted event from model 3

xifore.m1 <- NULL
len.chain <- length((m1.mcmc[,1])[[1]])
for (i in 1:len.chain){ 
  tem <- as.vector(m1.mcmc[i,][[1]])
  param <- tem[c(1,2,6)]
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  xifore.m1[i] <- rweibull(1,shape=alpha,scale=Z^(-1/alpha)*beta)
}

xifore.m2 <- NULL
len.chain <- length((m2.mcmc[,1])[[1]])
for (i in 1:len.chain){ 
  tem <- as.vector(m2.mcmc[i,][[1]])
  param <- tem[c(1,2,6)]
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  xifore.m2[i] <- rweibull(1,shape=Z*alpha,scale=Z^(1/Z/alpha)*beta^(1/Z))
}

xifore.m3 <- NULL
len.chain <- length((m3.mcmc[,1])[[1]])
for (i in 1:len.chain){ 
  tem <- as.vector(m3.mcmc[i,][[1]])
  param <- tem[c(1,2,10,7)]
  alpha <- param[1]
  beta <- param[2]
  Z <- param[3]
  Y <- param[4]
  xifore.m3[i] <- rweibull(1,shape=Z*alpha,scale=Z^(1/Z/alpha)*Y^(-1/Z/alpha)*beta^(1/Z))
}


postscript("LongRecV3AnalogWeib-M123-PredDistr.eps",paper="special",width=1.5*5*cos (35.4/180*pi)/0.612,
    height=0.5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,3),mar=c(4.5,4.5,1,1))
hist(log(xifore.m1),xlab="Waiting time (x1000 years)",xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     breaks=20,freq=F,xlim=c(log(0.0000001),log(max(xifore.m1))))
axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,1,10,100)),c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
     line=-0.5,cex.axis=1.5,lwd=0,lwd.ticks=1.5)

hist(log(xifore.m2),xlab="Waiting time (x1000 years)",xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     breaks=20,freq=F,xlim=c(log(0.0000001),log(max(xifore.m2))))
axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,1,10,100)),c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
     line=-0.5,cex.axis=1.5,lwd=0,lwd.ticks=1.5)

hist(log(xifore.m3),xlab="Waiting time (x1000 years)",xaxt="n",cex.lab=1.5,cex.axis=1.5,main="",
     breaks=20,freq=F,xlim=c(log(0.0000001),log(max(xifore.m1))))
axis(1,at=log(c(0.00001,0.0001,0.001,0.01,0.1,1,10,100)),c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
     line=-0.5,cex.axis=1.5,lwd=0,lwd.ticks=1.5)
dev.off()






gelman.diag(m1.mcmc)
gelman.diag(m2.mcmc)
gelman.diag(m3.mcmc)

options("scipen"=100, "digits"=4)
apply(cbind(p1c1,p1c2,p1c3),1,mean)
apply(cbind(p2c1,p2c2,p2c3),1,mean)
apply(cbind(p3c1,p3c2,p3c3),1,mean)
apply(cbind(pY3c1,pY3c2,pY3c3),1,mean)

quantile(xifore.m1,c(0.025,0.25,0.5,0.75,0.975))*1000
quantile(xifore.m2,c(0.025,0.25,0.5,0.75,0.975))*1000
quantile(xifore.m3,c(0.025,0.25,0.5,0.75,0.975))*1000



> quantile(xifore.m1,c(0.025,0.25,0.5,0.75,0.975))*1000
      2.5%        25%        50%        75%      97.5% 
    0.2895    57.8093   400.4483  1835.3919 17750.9801 
> quantile(xifore.m2,c(0.025,0.25,0.5,0.75,0.975))*1000
       2.5%         25%         50%         75%       97.5% 
    0.03943    51.52458   612.41658  3856.27271 52890.10436 
> quantile(xifore.m3,c(0.025,0.25,0.5,0.75,0.975))*1000
       2.5%         25%         50%         75%       97.5% 
    0.06963    37.74687   336.66937  1903.07247 26941.89257 



	
