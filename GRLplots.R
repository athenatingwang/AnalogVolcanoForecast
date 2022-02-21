# Quaternary data
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

col.tgr <- rep(1,length(tgr3$VEI))
col.tgr[tgr3$VEI==4.5] <- 2#grey(0.5)
col.ths <- rep(1,length(ths3$VEI))
col.ths[ths3$VEI==4.5] <- 2#grey(0.5)
col.pey <- rep(1,length(pey3$VEI))
col.pey[pey3$VEI==4.5] <- 2#grey(0.5)
pch.tgr <- rep(16,length(tgr3$VEI))
pch.tgr[tgr3$VEI==4.5] <- 17
pch.ths <- rep(16,length(ths3$VEI))
pch.ths[ths3$VEI==4.5] <- 17
pch.pey <- rep(16,length(pey3$VEI))
pch.pey[pey3$VEI==4.5] <- 17


## GVP Holocene data
source("ReadGVPdata.R")

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
  
Nind <- c(21,1:20,22:23)


## VEI versus time plot for all data
postscript("VEI-VS-time-geological-empirical.eps",paper="special",
    width=1.5*5*cos (35.4/180*pi)/0.612,
    height=1.5*5*7/4,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(7,4), mar=c(4.5, 4.5, 1, 1))

## VEI versus time plot for all data
i <- 1
for (kk in Nind){
  ttt <- get(paste('d',substr(a[kk],1,4),sep=""))
  colh <- rep("black",length(ttt$ti))
  colh[ttt$Evidence.Method..dating.=="Historical Observations"] <- grey(0.5)
  plot(ttt$ti,jitter(ttt$VEI,amount=0.3),xlab="time (year)",ylab="VEI",
       main=paste(aname[i],"(GVP)",sep=" "),pch=16,
       cex.lab=1.5,cex.axis=1.5,cex=0.8,col=colh,ylim=c(0,5),xlim=c(-10000,2018))
#  ttt <- get(paste('hisd',substr(a[kk],1,4),sep=""))
#  points(ttt$ti,jitter(ttt$VEI,amount=0.3),cex=0.8,pch=16,col=grey(0.5))
  i <- i+1
}

plot(1,1,type="n",axes=F)

plot(-tgr3$AgeyrBP,jitter(tgr3$VEI,amount=0.1),ylim=c(-0.1,5.1),
   main="Tongariro (geological)",xlab="time (x1000 year BP)",ylab="VEI",
   col=col.tgr,pch=pch.tgr,cex.lab=1.5,cex.axis=1.5,cex=1,xaxt="n")
axis(1,at=c(-15000,-10000,-5000,0),c(15,10,5,0),cex.axis=1.5)
text(-15000,5,"(a)",cex=1.8)

plot(-ths3$AgeyrBP,jitter(ths3$VEI,amount=0.1),ylim=c(-0.1,5.1),
   main="Three Sisters (geological)",xlab="time (x1000 year BP)",ylab="VEI",
   col=col.ths,pch=pch.ths,cex.lab=1.5,cex.axis=1.5,cex=1,xaxt="n")
axis(1,at=seq(-500000,0,100000),seq(500,0,-100),cex.axis=1.5)
text(-500000,5,"(b)",cex=1.8)

plot(-pey3$AgeyrBP,jitter(pey3$VEI,amount=0.1),ylim=c(-0.1,5.1),
   main="Puyehue-Cordon (geological)",xlab="time (x1000 year BP)",ylab="VEI",
   col=col.pey,pch=pch.pey,cex.lab=1.5,cex.axis=1.5,cex=1,xaxt="n")
axis(1,at=seq(-300000,0,100000),seq(300,0,-100),cex.axis=1.5)
text(-300000,5,"(c)",cex=1.8)
dev.off()





############### RESULTS ###############

#######################################
#### Quaternary record VEI 3+
#######################################
# R2jags
library(R2jags)
library(lattice)
source("denslines.R")
source("line2user.R")


# model 1
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm1-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm1mc-1010000i10000b.image"))
m1lv3.mcmc <- m1.mcmc

# model 2
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm2-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm2mc-1010000i10000b.image"))
m2lv3.mcmc <- m2.mcmc

# model 3
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm3-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI3+/TgriroLongRecV3Weibm3mc-1010000i10000b.image"))
m3lv3.mcmc <- m3.mcmc

## Probability that a VEI3+ eruption in the next 50 years
## 2070-1798=272
library(MCMCvis)
MCMCsummary(m3lv3.mcmc, 
          params = "xi",
          Rhat = TRUE,
          n.eff = TRUE,
          round = 2,
          func = function(x) ecdf(x)(272/1000),
          func_name = "ecdf50")


#######################################
#### Quaternary record VEI 4+
#######################################

# model 1
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm1-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm1mc-1010000i10000b.image"))
m1lv4.mcmc <- m1.mcmc

# model 2
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm2-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm2mc-1010000i10000b.image"))
m2lv4.mcmc <- m2.mcmc

# model 3
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm3-1010000i10000b.image"))
get(load("ResultsGRL/LongRecWeibsepZuse.jags-Final-VEI4+/TgriroLongRecWeibm3mc-1010000i10000b.image"))
m3lv4.mcmc <- m3.mcmc

## Probability that a VEI4+ eruption in the next 50 years
## 2500+2070-1950=2620
library(MCMCvis)
MCMCsummary(m3lv4.mcmc, 
          params = "xi",
          Rhat = TRUE,
          n.eff = TRUE,
          round = 2,
          func = function(x) ecdf(x)(2620/1000),
          func_name = "ecdf50")


#######################################
#### GVP record empirical analogs VEI 3+ eruptions, finding statistical analogs
#######################################

# model 1
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m1-110000i10000b.image"))
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m1mc-110000i10000b.image"))
m1gvpemp.mcmc <- m1.mcmc

# model 2		
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m2-110000i10000b.image"))
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m2mc-110000i10000b.image"))
m2gvpemp.mcmc <- m2.mcmc

# model 3
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m3-110000i10000b.image"))
get(load("ResultsGRL/GVPM0-M3-EmpAnalogV3/TgriroGVPEmpV3m3mc-110000i10000b.image"))
m3gvpemp.mcmc <- m3.mcmc


#######################################
#### GVP record statistical analogues VEI 3+
#######################################

# model 1
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m1-110000i10000b.image"))
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m1mc-110000i10000b.image"))
m1gvpstat.mcmc <- m1.mcmc

# model 2		
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m2-110000i10000b.image"))
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m2mc-110000i10000b.image"))
m2gvpstat.mcmc <- m2.mcmc

# model 3
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m3-110000i10000b.image"))
get(load("ResultsGRL/GVPM1-M3-StatsAnalogV3/TgriroGVPStatsV3m3mc-110000i10000b.image"))
m3gvpstat.mcmc <- m3.mcmc

## Probability that a VEI3+ eruption in the next 50 years
## 2070-1500=570
library(MCMCvis)
MCMCsummary(m3gvpstat.mcmc, 
          params = "xi",
          Rhat = TRUE,
          n.eff = TRUE,
          round = 2,
          func = function(x) ecdf(x)(570/1000),
          func_name = "ecdf50")


### Posteriors of Z and Y for all 3 models (Quaternary, empirical and Holocene analogues)
postscript("AnalogM123LongGVPpostdens.eps",paper="special",width=2*5*cos (35.4/180*pi)/0.612,
    height=0.5*2*5*3,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(6,4), mar=c(4.1, 5.1, 2, 1))

Mvol <- 23
temp <- rainbow(Mvol)
colsN <- temp 
colsN[1] <- 1
colsN[3] <- grey(0.4)


#### M1 GVP empirical analogues VEI 3+, finding statistical analogues
Mvol <- 23
m1.Zmcmc <- m1gvpemp.mcmc[,c(6,17,22:28,7:16,18:21)]
m1.sZmcmc <- m1gvpemp.mcmc[,4]

densplot(m1.Zmcmc[,1],xlim=c(0,5.5),ylim=c(0,3.2),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
  denslines(m1.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,3,"(a)",cex=1.5)

densplot(m1.sZmcmc,xlim=c(0,1.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3.65,"(b)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Empirical analogues (M1)', xpd=NA, 
     cex=1.8, font=2)
	
m2.Zmcmc <- m2gvpemp.mcmc[,c(6,17,22:28,7:16,18:21)]
m2.sZmcmc <- m2gvpemp.mcmc[,4]

densplot(m2.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m2.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,6.1,"(c)",cex=1.5)
densplot(m2.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(0.78,6.1,"(d)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Empirical analogues (M2)', 
     xpd=NA, cex=1.8, font=2)

m3.Ymcmc <- m3gvpemp.mcmc[,c(7,18,23:29,8:17,19:22)]
m3.Zmcmc <- m3gvpemp.mcmc[,c(30,41,46:52,31:40,42:45)]
m3.sYmcmc <- m3gvpemp.mcmc[,4]
m3.sZmcmc <- m3gvpemp.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5.5),ylim=c(0,4),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,3.7,"(e)",cex=1.5)

densplot(m3.sYmcmc,xlim=c(0,1.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3,"(f)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Empirical analogues (M3)', xpd=NA, 
     cex=1.8, font=2)

densplot(m3.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,5.7,"(g)",cex=1.5)
densplot(m3.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(0.77,5.5,"(h)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Empirical analogues (M3)', 
     xpd=NA, cex=1.8, font=2)


#### M1, M2, M3 statistical analogues with Holocene records for VEI 3+ eruptions
Mvol <- 18
m1.Zmcmc <- m1gvpstat.mcmc[,c(6,16:23,7:15)]
m1.sZmcmc <- m1gvpstat.mcmc[,4]

densplot(m1.Zmcmc[,1],xlim=c(0,5.5),ylim=c(0,2.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
  denslines(m1.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,2.32,"(i)",cex=1.5)

densplot(m1.sZmcmc,xlim=c(0,1.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3.1,"(j)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Holocene records (M1)', 
     xpd=NA, cex=1.8, font=2)


m2.Zmcmc <- m2gvpstat.mcmc[,c(6,16:23,7:15)]
m2.sZmcmc <- m2gvpstat.mcmc[,4]

densplot(m2.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m2.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,5.7,"(k)",cex=1.5)

densplot(m2.sZmcmc,xlim=c(0,0.8),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
for (i in 1:(ceiling(23/2)/2))
  legend(0.22,8-i*0.85,paste("i=",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.5,lwd=1.8)
for (i in (ceiling(23/2)/2+1):ceiling(23/2))
 legend(0.5,8-(i-ceiling(23/2)/2)*0.85,paste("i=",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.5,lwd=1.8)
text(0.77,7.4,"(l)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Holocene records (M2)', 
     xpd=NA, cex=1.8, font=2)


m3.Ymcmc <- m3gvpstat.mcmc[,c(7,17:24,8:16)]
m3.Zmcmc <- m3gvpstat.mcmc[,c(25,35:42,26:34)]
m3.sYmcmc <- m3gvpstat.mcmc[,4]
m3.sZmcmc <- m3gvpstat.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5.5),ylim=c(0,3),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,2.8,"(m)",cex=1.5)

densplot(m3.sYmcmc,xlim=c(0,1.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3,"(n)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Holocene records (M3)', 
     xpd=NA, cex=1.8, font=2)


densplot(m3.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,5.7,"(o)",cex=1.5)

densplot(m3.sZmcmc,xlim=c(0,0.8),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
for (i in (ceiling(23/2)+1):(ceiling(23/2)*2-ceiling(23/2)/2))
  legend(0.2,6.9-(i-ceiling(23/2))*0.8,paste("i=",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.5,lwd=1.8)
for (i in (ceiling(23/2)*2-ceiling(23/2)/2+1):23)
 legend(0.5,6.9-(i-(ceiling(23/2)*2-ceiling(23/2)/2)+1)*0.8,paste("i=",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.5,lwd=1.8)
text(0.77,6.7,"(p)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Holocene records (M3)', 
     xpd=NA, cex=1.8, font=2)



#statistical analogues with Quaternary records 

#### M1 Quaternary record VEI 3+
Mvol <- 3
m1.Zmcmc <- m1lv3.mcmc[,c(6:8)]
m1.sZmcmc <- m1lv3.mcmc[,4]

densplot(m1.Zmcmc[,1],xlim=c(0,5),ylim=c(0,2),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
  denslines(m1.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(4.8,1.9,"(q)",cex=1.5)

densplot(m1.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(14.8,0.95,"(r)",cex=1.5)
	
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Quaternary records (M1)', 
     xpd=NA, cex=1.8, font=2)


#### M2
m2.Zmcmc <- m2lv3.mcmc[,c(6:8)]
m2.sZmcmc <- m2lv3.mcmc[,4]

densplot(m2.Zmcmc[,1],xlim=c(0,2),ylim=c(0,4),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m2.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(1.9,3.7,"(s)",cex=1.5)

densplot(m2.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(6.6,1.9,"(t)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Quaternary records (M2)', 
     xpd=NA, cex=1.8, font=2)


#### M3
m3.Ymcmc <- m3lv3.mcmc[,c(7:9)]
m3.Zmcmc <- m3lv3.mcmc[,c(10:12)]
m3.sYmcmc <- m3lv3.mcmc[,4]
m3.sZmcmc <- m3lv3.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5),ylim=c(0,2),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(4.8,1.9,"(u)",cex=1.5)

densplot(m3.sYmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(9.6,1.06,"(v)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Quaternary records (M3)', 
     xpd=NA, cex=1.8, font=2)


densplot(m3.Zmcmc[,1],xlim=c(0,2),ylim=c(0,4.5),cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(1.9,4.2,"(w)",cex=1.5)

densplot(m3.sZmcmc,cex.axis=1.5,cex.lab=1.5,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(6.7,2.5,"(x)",cex=1.5)

text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=1, side=3), 'Statistical analogues Quaternary records (M3)', 
     xpd=NA, cex=1.8, font=2)

dev.off()



### Posteriors of Z and Y for model 3 (Quaternary, empirical and Holocene analogues)
####### Only Model 3 results
postscript("AnalogM3LongGVPpostdens.eps",paper="special",width=2*5*cos (35.4/180*pi)/0.612,
    height=0.5*5*3.5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(3,4), mar=c(4.5, 5.1, 2.3, 1))

Mvol <- 23
temp <- rainbow(Mvol)
colsN <- temp 
colsN[1] <- 1
colsN[3] <- grey(0.4)


#### M3 GVP empirical analogues VEI 3+, finding statistical analogues with Holocene records
Mvol <- 23

m3.Ymcmc <- m3gvpemp.mcmc[,c(7,18,23:29,8:17,19:22)]
m3.Zmcmc <- m3gvpemp.mcmc[,c(30,41,46:52,31:40,42:45)]
m3.sYmcmc <- m3gvpemp.mcmc[,4]
m3.sZmcmc <- m3gvpemp.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5.5),ylim=c(0,4),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,3.7,"(a)",cex=1.8)

densplot(m3.sYmcmc,xlim=c(0,1.5),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3,"(b)",cex=1.8)

densplot(m3.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,5.7,"(c)",cex=1.8)
densplot(m3.sZmcmc,cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(0.77,5.5,"(d)",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+25, side=2), 
     line2user(line=1, side=3), 'Empirical analogues (M3)', 
     xpd=NA, cex=1.8, font=2)


#### M3  statistical analogues with Holocene records for VEI 3+
Mvol <- 18
m3.Ymcmc <- m3gvpstat.mcmc[,c(7,17:24,8:16)]
m3.Zmcmc <- m3gvpstat.mcmc[,c(25,35:42,26:34)]
m3.sYmcmc <- m3gvpstat.mcmc[,4]
m3.sZmcmc <- m3gvpstat.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5.5),ylim=c(0,3),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(5.3,2.8,"(e)",cex=1.8)

densplot(m3.sYmcmc,xlim=c(0,1.5),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
text(1.46,3,"(f)",cex=1.8)

densplot(m3.Zmcmc[,1],xlim=c(0,3),ylim=c(0,6),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(2.8,5.7,"(g)",cex=1.8)

densplot(m3.sZmcmc,xlim=c(0,0.8),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
text(0.77,6.6,"(h)",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+25, side=2), 
     line2user(line=1, side=3), 'Statistical analogues with Holocene records (M3)', 
     xpd=NA, cex=1.8, font=2)



#### M3 Statistical analogues with Quaternary record VEI 3+
Mvol <- 3
#### M3
m3.Ymcmc <- m3lv3.mcmc[,c(7:9)]
m3.Zmcmc <- m3lv3.mcmc[,c(10:12)]
m3.sYmcmc <- m3lv3.mcmc[,4]
m3.sZmcmc <- m3lv3.mcmc[,5]

densplot(m3.Ymcmc[,1],xlim=c(0,5),ylim=c(0,2),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Y[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Ymcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(4.8,1.9,"(i)",cex=1.8)

densplot(m3.sYmcmc,cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Y]),ylab="Density",lwd=2)
for (i in 1:(ceiling(23/2)/2))
  legend(2,1.1-i*0.13,paste("i=",i,sep=""),bty="n",
         lty=i,col=colsN[i],cex=1.6,lwd=1.6)
for (i in (ceiling(23/2)/2+1):ceiling(23/2))
 legend(6,1.1-(i-ceiling(23/2)/2)*0.13,paste("i=",i,sep=""),
        bty="n",lty=i,col=colsN[i],cex=1.6,lwd=1.6)
text(9.6,1.06,"(j)",cex=1.8)

densplot(m3.Zmcmc[,1],xlim=c(0,2),ylim=c(0,4.5),cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(Z[i]),ylab="Density",lwd=2)
for (i in 2:Mvol)
denslines(m3.Zmcmc[,i],col=colsN[i],lty=i,lwd=2,show.obs = FALSE)
text(1.9,4.2,"(k)",cex=1.8)

densplot(m3.sZmcmc,cex.axis=1.8,cex.lab=1.8,
    main="",xlab=expression(sigma[Z]),ylab="Density",lwd=2)
for (i in (ceiling(23/2)+1):(ceiling(23/2)*2-ceiling(23/2)/2))
  legend(1,2.7-(i-ceiling(23/2))*0.32,paste("i=",i,sep=""),
         bty="n",lty=i,col=colsN[i],cex=1.6,lwd=1.6)
for (i in (ceiling(23/2)*2-ceiling(23/2)/2+1):23)
 legend(4,2.7-(i-(ceiling(23/2)*2-ceiling(23/2)/2)+1)*0.32,
        paste("i=",i,sep=""),bty="n",lty=i,col=colsN[i],cex=1.6,lwd=1.6)
text(6.7,2.5,"(l)",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+25, side=2), 
     line2user(line=1, side=3), 'Statistical analogues with Quaternary records (M3)', 
     xpd=NA, cex=1.8, font=2)

dev.off()




###############################################
####      Forecasting
###############################################
m1.maxfore <- max(as.matrix(m1gvpstat.mcmc[,5]),as.matrix(m1gvpemp.mcmc[,5]))
m2.maxfore <- max(as.matrix(m2gvpstat.mcmc[,5]),as.matrix(m2gvpemp.mcmc[,5]))
m3.maxfore <- max(as.matrix(m3gvpstat.mcmc[,6]),as.matrix(m3gvpemp.mcmc[,6]))


postscript("TongariroV3V4Analog-Forecast.eps",paper="special",
    width=1.8*5*cos (35.4/180*pi)/0.612,
    height=0.5*5*5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(4,3),mar=c(5.1,4.5,4.2,1))

#################################
## Forecast x[6] from statistical analogues with Holocene records for VEI 3+
m1.x14mcmc <- m1gvpstat.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(4,0.47,"(a)",cex=1.8)
text(3.4,0.34,"M1",cex=1.8)
text(3.4,0.26,"VEI 3+",cex=1.8)

m2.x14mcmc <- m2gvpstat.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5))+0.1,
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(2.8,0.85,"(b)",cex=1.8)
text(2.2,0.62,"M2",cex=1.8)
text(2.2,0.45,"VEI 3+",cex=1.8)

m3.x14mcmc <- m3gvpstat.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.995)),
     lty="dotted",lwd=5,col="red")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.995)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.995))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8,col.axis="red")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(3.8,0.46,"(c)",cex=1.8)
text(3.4,0.32,"M3",cex=1.8)
text(3.4,0.25,"VEI 3+",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+15, side=2), 
     line2user(line=3.2, side=3), 
     'Using statistical analogues with Holocene records (VEI 3+)', 
     xpd=NA, cex=1.8, font=2)


#################################
## Forecast x[6] from GVP empirical analogues VEI 3+
m1.x14mcmc <- m1gvpemp.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(3.8,0.46,"(d)",cex=1.8)
text(3.3,0.34,"M1",cex=1.8)
text(3.3,0.26,"VEI 3+",cex=1.8)

m2.x14mcmc <- m2gvpemp.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5))+0.1,
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(2.2,0.95,"(e)",cex=1.8)
text(1.85,0.7,"M2",cex=1.8)
text(1.85,0.52,"VEI 3+",cex=1.8)

m3.x14mcmc <- m3gvpemp.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),max(log(xx[,1,drop=TRUE]))),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.995)),
     lty="dotted",lwd=5,col="red")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     c("","",floor(quantile((xx[,1,drop=TRUE]),c(0.975))*1000+1500)),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.995)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.995))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8,col.axis="red")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025))-0.1,
     floor(quantile((xx[,1,drop=TRUE]),c(0.025))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5))+0.2,
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1500),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     c("",""),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25))-0.1,
     floor(quantile((xx[,1,drop=TRUE]),c(0.25))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75))+0.1,
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000+1500),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=0)
box()
text(6.9,0.38,"(f)",cex=1.8)
text(5.8,0.28,"M3",cex=1.8)
text(5.8,0.21,"VEI 3+",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+15, side=2), 
     line2user(line=3.2, side=3), 
     'Using empirical analogues (VEI 3+)', 
     xpd=NA, cex=1.8, font=2)



##################################
## Forecast x[14] from Statistical analogues with Quaternary records for VEI 3+
m1.x14mcmc <- m1lv3.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1798),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1798),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1798),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(6.7,0.28,"(g)",cex=1.8)
text(5.7,0.23,"M1",cex=1.8)
text(5.7,0.18,"VEI 3+",cex=1.8)

m2.x14mcmc <- m2lv3.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1798),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1798),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(6.7,0.215,"(h)",cex=1.8)
text(5.7,0.175,"M2",cex=1.8)
text(5.7,0.135,"VEI 3+",cex=1.8)

m3.x14mcmc <- m3lv3.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7),
     breaks=seq(log(0.22),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
#quantile((xx[,1,drop=TRUE]),c(0.05,0.25,0.5,0.75,0.95))*1000+1798
#  2066.034  2373.197  3257.082  6203.279 25217.844 
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000+1798),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000+1798),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000+1798),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(6.7,0.265,"(i)",cex=1.8)
text(5.7,0.215,"M3",cex=1.8)
text(5.7,0.165,"VEI 3+",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+15, side=2), 
     line2user(line=3.2, side=3), 
     'Using statistical analogues with Quaternary records (VEI 3+)', 
     xpd=NA, cex=1.8, font=2)



###################################
## Forecast x[9] from Statistical analogues with Quaternary records VEI 4+
m1.x14mcmc <- m1lv4.mcmc[,5]
xx <- as.matrix(m1.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),6.5),
     breaks=seq(min(log(xx)),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000-550),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000-550),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(6.2,0.52,"(j)",cex=1.8)
text(5.5,0.42,"M1",cex=1.8)
text(5.5,0.32,"VEI 4+",cex=1.8)

m2.x14mcmc <- m2lv4.mcmc[,5]
xx <- as.matrix(m2.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7.5),
     breaks=seq(min(log(xx)),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5)),
    floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000-550),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.75))*1000-550),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
box()
text(7.2,0.42,"(k)",cex=1.8)
text(6.4,0.34,"M2",cex=1.8)
text(6.4,0.26,"VEI 4+",cex=1.8)

m3.x14mcmc <- m3lv4.mcmc[,6]
xx <- as.matrix(m3.x14mcmc)
hist(log(xx[,1,drop=TRUE]),xlab="Forecast time of next eruption (year)",
     xaxt="n",cex.lab=1.8,cex.axis=1.8,main="",
     xlim=c(min(log(xx[,1,drop=TRUE])),7.5),
     breaks=seq(min(log(xx)),max(log(xx))+0.1,length.out=50),freq=F)
abline(v=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.25,0.5,0.75,0.975)),
     lty="dashed",lwd=2,col="blue")
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.025,0.5,0.975)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.025,0.5,0.975))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(1,at=quantile(log(xx[,1,drop=TRUE]),c(0.5))+0.2,
     floor(quantile((xx[,1,drop=TRUE]),c(0.5))*1000-550),
     line=-0,cex.axis=1.8,lwd=0,lwd.ticks=0)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.25,0.75)),
     floor(quantile((xx[,1,drop=TRUE]),c(0.25,0.75))*1000-550),
     line=-0.2,cex.axis=1.8,lwd=0,lwd.ticks=1.8)
axis(3,at=quantile(log(xx[,1,drop=TRUE]),c(0.75)),
     c(15553),
     line=-0.2,cex.axis=1.8,lwd=0,lw3d.ticks=1.8)
box()
text(7.2,0.48,"(l)",cex=1.8)
text(6.4,0.4,"M3",cex=1.8)
text(6.4,0.31,"VEI 4+",cex=1.8)

text(line2user(line=mean(par('mar')[c(2,4)])+15, side=2), 
     line2user(line=3.2, side=3), 
     'Using statistical analogues with Quaternary records (VEI 4+)', 
     xpd=NA, cex=1.8, font=2)

dev.off()












