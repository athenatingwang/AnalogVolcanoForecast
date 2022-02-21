## VEI 3+
tgr3 <- read.csv("QuaternaryRecords/Tongariro.csv",header=TRUE)
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

m0 = jags(data=data1, inits = vinits, parameters.to.save=params, n.chains = 3,
          n.iter = 110000,n.burnin=10000,n.thin=20,
          model.file="ModelM0Weib.jags")

m0.mcmc <- as.mcmc(m0)
summary(m0.mcmc)


m3.x14mcmc <- m0.mcmc[,4]
xx <- as.matrix(m3.x14mcmc)
floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975))*1000+1798)


save(m0,file="TgriroLongRecV3Weibm0-1010000i10000b.image")
save(m0.mcmc,file="TgriroLongRecV3Weibm0mc-1010000i10000b.image")

