source("ReadGVPdata.R")

## Using eruption records with minimum VEI 3

aname <- c("Tongariro","Asamayama","Avachinsky","Azumayama","Bezymianny",
"Calbuco","Chokaisan","Colima","Fujisan","Gamalama","Guagua Pichincha",
"Kuchinoerabujima","Miyakejima","Nevado del Ruiz","Pico de Orizaba",
"Popocatepetl","Puyehue-Cordon Caulle","Sheveluch","Suwanosejima",
"Tokachidake","Tolbachik","Tungurahua","Villarrica")

Nind <- c(21,1:20,22:23)

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
  ## tiV3: VEI 3+ eruptions only
  temp <- eval(parse(text=paste("tiV3",substr(a[kj],1,4),sep="")))
  len.volc[kj] <- length(temp)
  max.erp <- c(max.erp,temp[length(temp)])
  eval(parse(text=paste(substr(a[kj],1,4),'IE=diff(temp)/1000',sep="")))
  eval(parse(text=paste("y",substr(a[kj],1,4),
         '=(lastday[1]-temp[length(temp)])/1000',sep="")))
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
N1 <- c(1,N2[-length(N2)]+1) # The indices for the
                      # first data points from all the volcano in dat
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
m1gvpemp.mcmc <- m1.mcmc
summary(m1.mcmc)

save(m1,file="TgriroGVPEmpV3m1-110000i10000b.image")
save(m1.mcmc,file="TgriroGVPEmpV3m1mc-110000i10000b.image")

m1.Zmcmc <- m1.mcmc[,c(6,17,22:28,7:16,18:21)]
m1.sZmcmc <- m1.mcmc[,4]

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
  list("alpha" = 0.7, "beta" = 0.3,"sigmaZ" = 1, "Z"=rep(1,Mvol))
}

params <- c("alpha", "beta", "sigmaZ", "Z","xi[6]")

m2 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="ModelM2Weib.jags")

m2.mcmc <- as.mcmc(m2)
m2gvpemp.mcmc <- m2.mcmc
summary(m2.mcmc)

save(m2,file="TgriroGVPEmpV3m2-110000i10000b.image")
save(m2.mcmc,file="TgriroGVPEmpV3m2mc-110000i10000b.image")

m2.Zmcmc <- m2.mcmc[,c(6,17,22:28,7:16,18:21)]
m2.sZmcmc <- m2.mcmc[,4]

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

options("scipen"=100, "digits"=4)
apply(cbind(p2c1,p2c2,p2c3),1,mean)




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
m3gvpemp.mcmc <- m3.mcmc
summary(m3.mcmc)

save(m3,file="TgriroGVPEmpV3m3-110000i10000b.image")
save(m3.mcmc,file="TgriroGVPEmpV3m3mc-110000i10000b.image")

m3.x14mcmc <- m3.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.005,0.025,0.05,0.25,
      0.5,0.75,0.95,0.975,0.995))*1000+1500)

library("MCMCvis")
MCMCsummary(m3.mcmc,HPD=TRUE,hpd_prob=0.95,round=3)

m3.Ymcmc <- m3.mcmc[,c(7,18,23:29,8:17,19:22)]
m3.Zmcmc <- m3.mcmc[,c(30,41,46:52,31:40,42:45)]
m3.sYmcmc <- m3.mcmc[,4]
m3.sZmcmc <- m3.mcmc[,5]

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

