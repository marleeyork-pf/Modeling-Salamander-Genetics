# This script is for the implementation and analysis of past spatial and genomic
# cline models, and the new application of their joint modeling.
#Packages
#library(gghybrid);library(MASS)
library(misc3d);library(viridis)
library(viridisLite)
#library(hzar);library(HIest) 
library(optimx);library(scales)
library(png);library(RIdeogram)
library(hierfstat);library(pegas)
#library(arm);library(demodelr)
library(plotly);library(iterators)
library(foreach);library(parallel)
library(doParallel);library(plyr)

#Functions
source("./aux_funcs.R")
source("./marlee_aux.R")

#Sample data for monticola - 71 samples from Pyron et al. 2023
dat <- read.csv("./GBS/gbs_71_localities.csv",header=T,row.names=1)

#SNPs for monticola - 310 loci from Pyron et al. 2023, rescaled to monticola
snps <- read.csv("./GBS/monticola_310_SNP_freqs_rescaled.csv",header=T,row.names=1)

#Locus names - to be used later for locus-specific estimates
loci <- read.csv("./GBS/Locus_Names.csv",header=T,row.names=1)
head(loci)

#Match the loci to our SNPs
match(colnames(snps),rownames(loci))
loci$ID[match(colnames(snps),rownames(loci))]

#Put all the data together
dd.all <- data.frame(x=dat$dist,
                     x.s=trans.down(dat$dist),
                     x.g=dat$monticola,
                     n=2,
                     y=dat$monticola,
                     z=snps*2)
dd <- data.frame(na.omit(dd.all[,1:6]))#Just one SNP


##################################
#Plot a spatial cline - Admixture#
##################################

res.admixture <- multistart(parmat=pars.A,
                            fn=LL.admixture,
                            x=dd.all$x,n=dd.all$n,y=dd.all$y,
                            lower=c(min(dd.all$x)*1.1,0),
                            upper=c(max(dd.all$x)*1.1,max(dd$x)*2.2),
                            method="L-BFGS-B",
                            control=list(maximize=TRUE))

est.admixture <- res.admixture[which.max(res.admixture$value),]
est.admixture$aic <- -est.admixture$value*2+4

#plot
plot(x=dd.all$x,y=dd.all$y,type="n",xlab=expression(paste(italic(S))),
     ylab=expression(paste(italic(p[i])," beta")))
curve(Barton.p.noTails(x,unlist(est.admixture[,1:2])),from=min(dd.all$x),to=max(dd.all$x),col="red",lwd=2,add=T)
curve(Barton.p.noTails(x,c(0,0)),from=min(dd.all$x),to=max(dd.all$x),col="darkgrey",lty=2,lwd=2,add=T)
points(dd.all[,c("x","y")],pch=21,bg=c("#1f77b4","#ff7f0e")[round(dd.all$y)+1],cex=2)
est.admixture


#####################################
#Plot a joint cline - Heterozygosity#
#####################################
dens <- kde3d(dd$x, dd$y, -dd$z/2)  #overrode default bandwidth
filled.contour(dens$x, dens$y, dens$d[,,3], color = function(n) plasma(n),
               xlim=c(-500,600),ylim=c(-0.05,1.05),
               plot.title = title(main = "Heterozygosity Surface (spatial/genomic)",
                                  xlab = "Cline Distance", ylab = "Genomic Ancestry"),
               plot.axes = {axis(1, seq(-800, 800, by = 100))
                 axis(2, seq(0, 1, by = 0.1)) 
                 points(dd[,c("x","y")],pch=21,bg=rev(viridis(3))[factor(-dd$z/2)],cex=2)},
               key.title = title(main = "pA"),
               key.axes = axis(4, seq(0, 1, by = 0.1)))
legend(200,0.4,legend=c("AA","Aa","aa"),fill=rev(viridis(3)),cex=2)


##########################
###Fit 2D spatial clines #
##Szymura and Barton 1986#
##Fit with/without tails #
##########################

###EXAMPLE FOR RESIDUALS

#Get cline with no tails
test.noTails <- fit.noTails(dd$x.s,dd$n,dd$z/2)

#Plot
plot(x=dd$x,y=dd$z/2,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," beta")),main="Barton Spatial Cline")
curve(Barton.p.noTails(x,c(0,0)),from=min(dd$x),to=max(dd$x),col="darkgrey",lty=2,lwd=2,add=T) #NULL
curve(Barton.p.noTails(x,c(trans.up(test.noTails$c,dd$x),test.noTails$w*diff(range(dd$x)))),
      from=min(dd$x),to=max(dd$x),col="black",lwd=5,add=T) #NULL
points(dd$x,dd$z/2,pch=21,bg=viridis(3)[dd$z+1],cex=2)

#Get the probability densities
dd$p <- Barton.p.noTails(dd$x,c(trans.up(test.noTails$c,dd$x),test.noTails$w*diff(range(dd$x))))

#Sort everything
dd <- dd[order(dd$x),]

#Get the distances
x <- dd$x

#Get the genotypes
y <- dd$z

#Get the probabilities
p <- dd$p


##############################################
### DETERMINING RESIDUAL ANALYSIS APPROACH ###
##############################################

# ACCURACY (same idea as rounding and calculating error)
# Plot of discrete residuals.
y.expected <- round(p*2)
plot(x=x, y=y.expected, type='l', main="Discrete Residuals: Accuracy", xlab="S", ylab="Allele Count", 
     col="hotpink1", pch=20, font.main=4, lwd=3, ylim=c(-.1, 2.1))
points(x=x, y=y, col="darkmagenta")
for (i in 1:length(x)) {
  segments(x0=x[i], y0=y[i], x1=x[i], y1=y.expected[i], col="lemonchiffon3")
}

legend(125,.5,legend=c("Expected (rounded)", "Observed"),fill=c("hotpink1","darkmagenta"))

accuracy <- sum(y.expected==y) / length(y)

# RAW RESIDUALS
# plot of raw residuals, binned plot, and calculating MSE
plot(x=x, y=y/2, main="Continuous Residuals", 
     xlab="S", ylab="Allele Frequency", col="hotpink1",
     pch=20, font.main=4, lwd=3)
curve(Barton.p.noTails(x,c(trans.up(test.noTails$c,x),test.noTails$w*diff(range(x)))),from=min(x),to=max(x),col="darkmagenta",lwd=3,add=T)
for (i in 1:length(x)) {
  segments(x0=x[i], y0=y[i]/2, x1=x[i], y1=p[i], col="lemonchiffon3")
}
legend(125,.4,legend=c("Observed", "Expected"),fill=c("hotpink1","darkmagenta"))


mse <- sum(((y/2)-p)^2)

# DEVIANCE RESIDUALS
point.dev <- sqrt(2*(-log(1-abs((y/2)-p))))
for (i in 1:length(y)){
  if (y[i]/2 < p[i]) {
    point.dev[i] = -point.dev[i]
  }
}
residual.deviance <- sum(point.dev^2)

# plot of deviance residuals
plot(x=x, y=point.dev, main="Deviance Residuals",
     xlab="S", ylab="Residual", col="black", font.main=4)
abline(h=0,col="red")

# Error measures calculated so far
measures <- data.frame(c("Accuracy", "MSE", "Res Deviance"), c(accuracy, mse, residual.deviance))
colnames(measures) <- c("Measure", "Value")
measures

# SPLITTING COUNTS FOR RESIDUAL DEVIANCE
allele1 <- rep(0,length(y))
allele2 <- rep(0,length(y))
for (i in 1:length(y)){
  if (y[i]==1){
    allele1[i]=1
    allele2[i]=0
  }
  if (y[i]==2){
    allele1[i]=1
    allele2[i]=1
  }
}
alleles <- data.frame(allele1, allele2)
alleles

# plots of residuals and deviance residuals separately

par(mfrow=c(1,2))

plot(x=x,y=allele1, xlab="S", ylab="Allele 1", main="1st Allele Presence")
curve(Barton.p.noTails(x,c(trans.up(test.noTails$c,x),test.noTails$w*diff(range(x)))),from=min(x),to=max(x),col="darkmagenta",lwd=3,add=T)
plot(x=x,y=allele2, xlab="S", ylab="Allele 2", main="2nd Allele Presence")
curve(Barton.p.noTails(x,c(trans.up(test.noTails$c,x),test.noTails$w*diff(range(x)))),from=min(x),to=max(x),col="darkmagenta",lwd=3,add=T)

#Calculating point residual seperately, then summing and squaring to get residual deviance.
point.dev1 <- sqrt(2*(-log(1-abs(allele1-p))))
for (i in 1:length(y)){
  if (allele1[i] < p[i]) {
    point.dev1[i] = -point.dev1[i]
  }
}

point.dev2 <- sqrt(2*(-log(1-abs(allele2-p))))
for (i in 1:length(y)){
  if (allele2[i] < p[i]) {
    point.dev2[i] = -point.dev2[i]
  }
}
resdev <- sum((point.dev1+point.dev2)^2)
resdev

# PLOT COMPARING RESIDUAL DEVIANCE WHEN SUMMED AND NOT SUMMED
par(mfrow=c(1,2))
plot(x=x, y=point.dev, main="Deviance Residuals",
     xlab="S", ylab="Residual", font.main=4,
     col=ifelse(y==1, "red", "black"))
abline(h=0, col="black")
plot(x=x, y=point.dev1+point.dev2, col=ifelse(y==1, "red", "black"), 
     main="Deviances (Summed Prior)", font.main=4, xlab="S", ylab="Residuals (Calculated Seperate)")
abline(h=0, col="black")

###TAILS & FREQUENCIES
# 
# let's do two things
# first, let's implement the likelihood for a genomic cline
# with continuous frequencies on [0,1] for the response variable
# this will be like what you already did for the admixture cline
# but for the genomic cline Fitzpatrick function

# here's an example

#Likelihood for continuous frequencies
likelihood.C <- function(p,y,n){((1-y)*log((1-p)/(1-y))+y*log(p/y))*n}

#Likelihood for discrete genotypes
likelihood.D <- function(p,y,n){dbinom(y,n,p,log=TRUE)}

#New likelihood for continuous
likelihood.C.freq <- function(pars.G,x,y,p)
{log((x^pars.G[2]/((x^pars.G[2])+(1-x)^pars.G[2]*pars.G[1]^exp(1)))*(p)+
       (1-(x^pars.G[2]/(x^pars.G[2]+(1-x)^pars.G[2]*pars.G[1]^exp(1))))*(1-p))}

#Fitzpatrick logit-logistic genomic cline 
pars.G <- data.frame(mu=runif(100,-100,100),
                     nu=runif(100,0,1000))

#PDF for the Fitzpatrick logit-logistic genomic cline
Fitz.p <- function(x,pars.G){x^pars.G[2]/(x^pars.G[2]+(1-x)^pars.G[2]*exp(pars.G[1]))}

#Fitzpatrick logit-logistic likelihood -- New likelihood implemented
LL.Fitz <- function(pars.G,x,y)
{p <- Fitz.p(x,pars.G)
p <- replace(p, p <= 0, .Machine$double.xmin)
p <- replace(p, p >= 1, 1 - .Machine$double.neg.eps)
LL <- likelihood.C.freq(pars.G,x,y,p)
return(log(sum(LL)))}


#Fitzpatrick logit-logistic fit
fitFitz <- function(pars.G,x,n,y)
{res.Fitz <- multistart(par=pars.G,
                        fn=LL.Fitz,
                        x=x,y=y,n=n,
                        lower=c(-Inf,0.0),
                        upper=c(Inf,Inf),
                        method="L-BFGS-B",
                        control=list(maximize=TRUE))
est.Fitz <- round(res.Fitz[which.max(res.Fitz$value),],2)
est.Fitz}

#Estimate Fitzpatrick genomic cline
est.Fitz <- fitFitz(pars.G,dd$x.g,dd$n,dd$z)
est.Fitz

#plot
plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," Fitz")))
curve(Fitz.p(x,c(0,1)),from=0,to=1,col="darkgrey",lty=2,lwd=2,add=T)
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,col="red",lwd=2,add=T)
points(dd$x.g,dd$z/2,pch=21,bg=viridis(3)[dd$z+1],cex=2)

#the likelihood for this is given on the bottom of page 5 of Bailey 2023
#for some reason, he does not have the equations numbered!
#I think you can drop the prior term for ln(pi0) given at the end

dd <- data.frame(na.omit(dd.all[,c(1:5,315)]))

#PDF for the Barton spatial cline, right tail
Barton.p.rTails <- function(x,pars.R)
{ifelse(pars.R[3] < x - pars.R[1], 1 - 1/(1 + 
                                            exp(4 * pars.R[3]/pars.R[2])) * exp(-(pars.R[4]/(1 + exp(-4 * pars.R[3]/pars.R[2]))) * 
                                                                                  ((x - pars.R[1]) * 4/pars.R[2] - 4 * pars.R[3]/pars.R[2])), 
        1/(1 + exp(-((x - pars.R[1]) * 4/pars.R[2]))))}


#PDF for the Barton spatial cline, left tail
Barton.p.lTails <- function(x,pars.L)
{ifelse(pars.L[3] < -(x - pars.L[1]), 1/(1 + 
                                           exp(4 * pars.L[3]/pars.L[2])) * exp(pars.L[4]/(1 + exp(-4 * pars.L[3]/pars.L[2])) * 
                                                                                 ((x - pars.L[1]) * 4/pars.L[2] + 4 * pars.L[3]/pars.L[2])), 
        1/(1 + exp(-((x - pars.L[1]) * 4/pars.L[2]))))}


# Parameters for Fitz with right Barton tail
pars.Fitz.R <- data.frame(mu=runif(100,-100,100),
                          nu=runif(100,0,1000),
                          d=runif(100,0,2),
                          t=runif(100))

pars.Fitz.L <- data.frame(mu=runif(100,-100,100),
                          nu=runif(100,0,1000),
                          d=runif(100,0,2),
                          t=runif(100))

#PDF for Fitz, right tail
Fitz.p.rTails <- function(x,pars.Fitz.R)
{ifelse(pars.Fitz.R[3] < x - pars.Fitz.R[1], 1 - 1/(1 + 
                                                      exp(4 * pars.Fitz.R[3]/pars.Fitz.R[2])) * exp(-(pars.Fitz.R[4]/(1 + exp(-4 * pars.Fitz.R[3]/pars.Fitz.R[2]))) * 
                                                                                                      ((x - pars.Fitz.R[1]) * 4/pars.Fitz.R[2] - 4 * pars.Fitz.R[3]/pars.Fitz.R[2])),
        x^pars.Fitz.R[2]/(x^pars.Fitz.R[2]+(1-x)^pars.Fitz.R[2]*exp(pars.Fitz.R[1])))}

#PDF for Fitz, left tail
Fitz.p.lTails <- function(x,pars.Fitz.L)
{ifelse(pars.Fitz.L[3] < -(x - pars.Fitz.L[1]), 1/(1 + 
                                                     exp(4 * pars.Fitz.L[3]/pars.Fitz.L[2])) * exp(pars.Fitz.L[4]/(1 + exp(-4 * pars.Fitz.L[3]/pars.Fitz.L[2])) * 
                                                                                                     ((x - pars.Fitz.L[1]) * 4/pars.Fitz.L[2] + 4 * pars.Fitz.L[3]/pars.Fitz.L[2])), 
        x^pars.Fitz.L[2]/(x^pars.Fitz.L[2]+(1-x)^pars.Fitz.L[2]*exp(pars.Fitz.L[1])))}

#Likelihood function for Fitz, right tail
LL.Fitz.rTails <- function(pars.Fitz.R,x,y)
{p <- Fitz.p.rTails(x,pars.Fitz.R)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C.freq(pars.Fitz.R,x,y,p)
return(sum(LL))}

#Likelihood function for Fitz, left tail
LL.Fitz.lTails <- function(pars.Fitz.L,x,y)
{p <- Fitz.p.lTails(x,pars.Fitz.L)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C.freq(pars.Fitz.L,x,y,p)
return(sum(LL))}

#Fitting Fitz, right tail
fit.Fitz.rTails <- function(pars.Fitz.R,x,y)
{res.rTails <- multistart(par=pars.Fitz.R,
                          fn=LL.Fitz.rTails,
                          x=x,y=y,
                          lower=c(-Inf,0,0,0),
                          upper=c(Inf,Inf,2.2,1),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.rTails <- round(res.rTails[which.max(res.rTails$value),],2)
return(est.rTails)}

#Fitting Fitz, left tail
fit.Fitz.lTails <- function(pars.Fitz.L,x,y)
{res.lTails <- multistart(par=pars.Fitz.L,
                          fn=LL.Fitz.lTails,
                          x=x,y=y,
                          lower=c(-Inf,0,0,0),
                          upper=c(Inf,Inf,2.2,1),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.lTails <- round(res.lTails[which.max(res.lTails$value),],2)
est.lTails}

# Fitting the Fitzpatrick function
est.Fitz <- fitFitz(pars.G,dd$x.g,dd$n,dd$z)
est.Fitz.rTails <- fit.Fitz.rTails(pars.Fitz.R,dd$x.g,dd$z)
est.Fitz.lTails <- fit.Fitz.lTails(pars.Fitz.L,dd$x.g,dd$z)

est.Fitz
est.Fitz.rTails
est.Fitz.lTails

par(mfrow=c(1,1))
plot(x=dd$x.g,y=dd$z/2, main="Fitzpatrick with Tails", xlab="S", ylab="pA")
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,col="blue",add=TRUE, lwd=3)
curve(Fitz.p.rTails(x,unlist(est.Fitz.rTails[,1:4])),from=0,to=1,col="red",add=TRUE, lwd=3)
curve(Fitz.p.lTails(x,unlist(est.Fitz.lTails[,1:4])),from=0,to=1,col="orange",add=TRUE,lwd=3)
legend(.7,.4,legend=c("No tails","Right Tail" ,"Left tail"),fill=c("blue","red","orange"))


##########################
#Joint Clines for Ringers#
##########################
source("./aux_funcs.R")#Reset your functions!

dd <- data.frame(na.omit(dd.all[,c(1:5,19)]))#Just one SNP

#Why is this returning wild results? It seems like it should fix
#center and width at 0,0. Can you try other optimization routines in R?
#Or maybe MCMC or something? Or is the likelihood surface messed up?

#Spatial cline
test.noTails <- fit.noTails(dd$x.s,dd$n,dd$z)
test.noTails$c <- trans.up(test.noTails$c,dd.all$x)
test.noTails$w <- test.noTails$w*diff(range(dd.all$x))

# Plot
plot(x=dd$x,y=dd$z,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," beta")))
curve(Barton.p.noTails(x,c(0,0)),from=min(dd$x),to=max(dd$x),col="darkgrey",lty=2,lwd=2,add=T) #NULL
curve(Barton.p.noTails(x,unlist(test.noTails[1:2])),from=min(dd$x),to=max(dd$x),col="black",lwd=5,add=T) #NULL
points(dd$x,dd$z,pch=21,bg=viridis(3)[dd$z+1],cex=2)


# Genomic cline
est.Fitz <- fitFitz(pars.G,dd$x.g,dd$n,dd$z)

# Plot
plot(0:1,0:1,type="n",xlab=expression(paste(italic(S))),ylab=expression(paste(italic(p[i])," Fitz")))
curve(Fitz.p(x,c(0,1)),from=0,to=1,col="darkgrey",lty=2,lwd=2,add=T)
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,col="red",lwd=2,add=T)
points(dd$x.g,dd$z,pch=21,bg=viridis(3)[dd$z+1],cex=2)


# Joint Spatial/Genomic Cline, No Tails#
# Estimate joint
est.J <- fitJoint.noTails(pars.J,dd$x.s,dd$x.g,dd$n,dd$z)
est.J

# Plot spatial
plot(x=dd$x,y=dd$z,type="n",xlab=expression(paste(italic(S))),
     ylab=expression(paste(italic(p[i])," beta")))
curve(Barton.p.noTails(x,c(0,0)),from=min(dd$x),to=max(dd$x),col="darkgrey",lty=2,lwd=2,add=T)
curve(Barton.p.noTails(x,unlist(test.noTails[1:2])),from=min(dd$x),to=max(dd$x),col="red",lwd=2,add=T) #NULL
curve(Barton.p.noTails(x,unlist(est.J[,1:2])),from=min(dd$x),to=max(dd$x),col="blue",lwd=2,add=T)
points(dd$x,dd$z,pch=21,bg=viridis(3)[dd$z+1],cex=2)

#plot genomic
plot(x=0:1,y=0:1,type="n",xlab=expression(paste(italic(S))),
     ylab=expression(paste(italic(p[i])," Fitz")))
curve(Fitz.p(x,c(0,1)),from=0,to=1,col="darkgrey",lty=2,lwd=2,add=T)
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,col="red",lwd=2,add=T)
curve(Fitz.p(x,unlist(est.J[,3:4])),from=0,to=1,col="blue",lwd=2,add=T)
points(dd$x.g,dd$z,pch=21,bg=viridis(3)[dd$z+1],cex=2)

#######################################################
### Implementing continuous likelihood into ringers ###
#######################################################
# Resetting functions
source("./aux_funcs.R")

# Piecewise continous likelihood function
likelihood.C <- function(p,y,n){
  df <- data.frame(p,y)
  df$LL <- with(df, ifelse(df$y==0,log(1-p)*n,
                           ifelse(df$y==1,log(p)*n,
                                  ((1-y)*log((1-p)/(1-y))+y*log(p/y))*n)))
  return(df$LL)
}

dd <- data.frame(na.omit(dd.all[,c(1:5,9)]))#Just one SNP

# Estimating/plotting continuous genomic
est.Fitz.cont <- fitFitz.cont(pars.G,dd$x.g,dd$z,dd$n)
est.Fitz <- fitFitz(pars.G,dd$x.g,dd$n,dd$z)

curve(Fitz.p(x,unlist(est.Fitz.cont[,1:2])),from=0,to=1,col="blue",lwd=2,main="Fitz")
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,col="red",lwd=2,add=T)
points(dd$x.g,dd$z,pch=21,bg=viridis(3)[factor(dd$z)],cex=2)
legend(.7,.4,legend=c("Continuous", "Discrete"),fill=c("blue","red"))

# Estimating/plotting spatial with piecewise continuous
est.noTails.cont <- fit.noTails.cont(dd$x.s,2,dd$z)
est.noTails.cont$c <- trans.up(est.noTails.cont$c,dd.all$x)
est.noTails.cont$w <- est.noTails.cont$w*diff(range(dd.all$x))

test.noTails <- fit.noTails(dd$x.s,dd$n,dd$z)
test.noTails$c <- trans.up(test.noTails$c,dd.all$x)
test.noTails$w <- test.noTails$w*diff(range(dd.all$x))

plot(dd$x,dd$z,pch=21,bg=viridis(3)[factor(dd$z)],cex=2,xlab="x",ylab="pA")
curve(Barton.p.noTails(x,unlist(est.noTails.cont[,1:2])),from=min(dd$x),to=max(dd$x),col="blue",lwd=2,add=T)
curve(Barton.p.noTails(x,unlist(test.noTails[1:2])),from=min(dd$x),to=max(dd$x),col="red",lwd=2,add=T)
legend(200,.3,legend=c("Continuous", "Discrete"),fill=c("blue","red"))

# estimating/plotting Joint
est.J.cont <- fitJoint.noTails.cont(pars.J,dd$x.s,dd$x.g,dd$n,dd$z)
est.J <- fitJoint.noTails(pars.J,dd$x.s,dd$x.g,dd$n,dd$z)

# genomic joint
curve(Fitz.p(x,unlist(est.J[,3:4])),from=0,to=1,col="red",lwd=2,ylim=c(-.1,1.1),main="Joint Fitz")
curve(Fitz.p(x,unlist(est.J.cont[,3:4])),from=0,to=1,col="blue",lwd=2,add=T)
curve(Fitz.p(x,unlist(est.Fitz.cont[,1:2])),from=0,to=1,col="green",lwd=2,add=T)
points(dd$x.g,dd$z,pch=21,bg=viridis(3)[factor(dd$z)],cex=2)
legend(.7,.4,legend=c("Continuous Joint", "Discrete Joint", "Fitz"),fill=c("blue","red","green"))

# spatial joint
curve(Barton.p.noTails(x,unlist(est.J[,1:2])),from=min(dd$x),to=max(dd$x),col="red",lwd=2,add=F,ylim=c(-.1,1.1),main="Joint Barton")
curve(Barton.p.noTails(x,unlist(est.J.cont[,1:2])),from=min(dd$x),to=max(dd$x),col="blue",lwd=2,ylim=c(-.1,1.1),add=T)
curve(Barton.p.noTails(x,unlist(est.noTails.cont[,1:2])),from=min(dd$x),to=max(dd$x),col="green",lwd=2,add=T)
points(dd$x,dd$z,pch=21,bg=viridis(3)[factor(dd$z)],cex=2)
legend(200,.2,legend=c("Continuous Joint", "Discrete Joint","Barton"),fill=c("blue","red","green"))

# stepping back for spatial
curve(Barton.p.noTails(x,unlist(est.J[,1:2])),from=-1000,to=1000,col="red",lwd=2,add=F,ylim=c(-.1,1.1))
curve(Barton.p.noTails(x,unlist(est.J.cont[,1:2])),from=-1000,to=1000,col="blue",lwd=2,ylim=c(-.1,1.1),add=T)
curve(Barton.p.noTails(x,unlist(est.noTails.cont[,1:2])),from=-1000,to=1000,col="green",lwd=2,add=T)
points(dd$x,dd$z,pch=21,bg=viridis(3)[factor(dd$z)],cex=2)
legend(200,.2,legend=c("Continuous Joint", "Discrete Joint","Barton"),fill=c("blue","red","green"))

##########################
# PDF's FOR BARTON PARAM #
##########################
likelihood.C <- function(p,y,n){
  df <- data.frame(p,y)
  df$LL <- with(df, ifelse(df$y==0,log(1-p)*n,
                           ifelse(df$y==1,log(p)*n,
                                  ((1-y)*log((1-p)/(1-y))+y*log(p/y))*n)))
  return(df$LL)
}


Barton.w.posterior <- function(c,x.s,y){
  
  w <- seq(0.000001,2,.02)
  LL.fun <- function(w){
    p <- Barton.p.noTails(x.s,c(c,w))
    p <- replace(p,p<=0,.Machine$double.xmin)
    p <- replace(p,p>=1,1-.Machine$double.neg.eps)
    return(sum(mapply(FUN=likelihood.C,p,y,2)))
  }
  LL <- mapply(FUN=LL.fun,w)
  #posterior <- LL / sum(LL)
  #return(posterior)
  return(-LL)
}

Barton.c.posterior <- function(w,x.s,y){
  c <- seq(0.000001,1,.01)
  LL.fun <- function(c){
    p <- Barton.p.noTails(x.s,c(c,w))
    p <- replace(p,p<=0,.Machine$double.xmin)
    p <- replace(p,p>=1,1-.Machine$double.neg.eps)
    return(sum(mapply(FUN=likelihood.C,p,y,2)))
  }
  LL <- mapply(FUN=LL.fun,c)
  #posterior <- LL / sum(LL)
  return(-LL)
}

dd <- data.frame(na.omit(dd.all[,c(1:5,6)]))

# plotting the spatial cline
par(mfrow=c(1,1))
est.noTails <- fit.noTails(dd$x.s,2,dd$z)
est.noTails$c <- trans.up(est.noTails$c,dd$x)
est.noTails$w <- est.noTails$w*diff(range(dd$x))
curve(Barton.p.noTails(x,unlist(est.noTails[,1:2])),from=min(dd$x),to=max(dd$x),col="green",lwd=2,ylim=c(-.1,1.1))
points(dd$x,dd$z,lwd=7,col="black")

# plotting the posterior distributions of w and c
# still not getting laplacian distribution
par(mfrow=c(1,2))
plot(w,unlist(Barton.w.posterior(.17,dd$x.s,dd$z)),main="Width Posterior",xlab="w",ylab="Posterior",pch=1)
plot(c,unlist(Barton.c.posterior(.84,dd$x.s,dd$z)),main="Center Posterior",xlab="c",ylab="Posterior",pch=1)

# Likelihood surface

# METHOD 1: multiplying posteriors of c and w
# make matrix of the values
c <- seq(.000001,1,.01)
w <- seq(0.000001,2,.02)

Matrix = matrix(0, nrow = 100, ncol = 100)
LLw <- as.data.frame(Matrix)
LLc <- as.data.frame(Matrix)

for (i in 1:length(c)){
  LLw[,i] <- unlist(Barton.w.posterior(c[i],dd$x.s,dd$z))
}

for (i in 1:length(w)){
  LLc[i,] <- unlist(Barton.c.posterior(w[i],dd$x.s,dd$z))
}

LL <- LLw + LLc
plot_ly(x=w,y=c,z=unlist(LL),type="contour")

# METHOD 2
# trying to solve for LL(W and C) = LL(W|C)*LL(C)
# calculating LL(C) = LL(C|W1) + ... + LL(C|W2)
c <- seq(.000001,1,.01)
w <- seq(0.000001,2,.02)

P.C <- 0
for (i in 1:length(w)){
  P.C <- P.C + Barton.c.posterior(w[i],dd$x.s,dd$z)
}

# calculating LL = LL(W|C)*LL(C)
Matrix = matrix(0, nrow = 100, ncol = 100)
LL <- as.data.frame(Matrix)
for (i in 1:length(c)){
  P.WC <- Barton.w.posterior(c[i],dd$x.s,dd$z)
  P.WandC <- P.WC+P.C
  LL[,i] <- unlist(P.WandC)
}

plot_ly(x=c,y=w,z=unlist(LL),type="contour")

# METHOD 3: directly calculating the joint likelihood?
# Results of this are equivalent to the previous method

c <- seq(.000001,1,.01)
w <- seq(0.000001,2,.02)
pars <- data.frame(c,w)

Matrix = matrix(0, nrow = 100)
LL <- as.data.frame(Matrix)


for (i in 1:length(c)){
  for (j in 1:length(w)){
    p <- Barton.p.noTails(dd$x.s,unlist(c(pars[i,1],pars[j,2])))
    p <- replace(p,p<=0,.Machine$double.xmin)
    p <- replace(p,p>=1,1-.Machine$double.neg.eps)
    
    plik <- 0
    for (k in 0:length(p)){
      plik[k] <- likelihood.C(p[k],dd$z[k],2)
    }
    LL[i,j] <- sum(plik)
  }
}

plot_ly(x=c,y=w,z=unlist(-LL),type="contour")

# Implementing pmin/pmax into the spatial cline

Barton.p.asy <- function(x,z,pars){
  min(z) + (max(z) - min(z)) * Barton.p.noTails(x,pars)
}

LL.noTails.asy <- function(pars.N,x,n,y,z)
{p <- Barton.p.asy(x,y,pars.N)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,z,n)
return(sum(LL))}

fit.noTails.asy <- function(x.s,n,y,z)
{res.noTails <- multistart(parmat=pars.N,
                           fn=LL.noTails.asy,
                           x=x.s,n=n,y=y,z=z,
                           lower=c(-0.1,0),
                           upper=c(1.1,2.2),
                           method="L-BFGS-B",
                           control=list(maximize=TRUE))
est.noTails <- round(res.noTails[which.max(res.noTails$value),],2)
est.noTails}
dd <- data.frame(na.omit(dd.all[,c(1:5,14)]))

# Trying on a normal snps - comes out identical

est.noTails <- fit.noTails(dd$x.s,dd$n,dd$z)
est.noTails.asy <- fit.noTails.asy(dd$x.s,dd$n,dd$y,dd$z)

plot(dd$x.s,dd$z,main="Asymtoping Barton Model",xlab="x",ylab="pA")
curve(Barton.p.noTails(x,unlist(est.noTails[,1:2])),from=0,to=1,add=TRUE,col="blue")
curve(Barton.p.asy(x,dd$z,unlist(est.noTails[,1:2])),from=0,to=1,add=TRUE,col="red") # identical

# Trying in the admixture cline
est.admixture <- fit.noTails(dd.all$x.s,dd.all$n,dd.all$y)
est.admixture.asy <- fit.noTails.asy(dd.all$x.s,dd.all$n,dd.all$y,dd.all$y)

plot(dd.all$x.s,dd.all$y,main="Asymptoping Barton Model",xlab="x",ylab="pA")
curve(Barton.p.noTails(x,unlist(est.admixture[,1:2])),from=0,to=1,add=TRUE,col="blue")
curve(Barton.p.asy(x,dd$y,unlist(est.admixture.asy[,1:2])),from=0,to=1,add=TRUE,col="red") # identical

###################################
####### D. Muller Testing #########
###################################
# Resetting functions and data
source("./aux_funcs.R")
samples <- read.csv("./GBS/gbs_71_localities.csv",header=T,row.names=1)
xy <- samples[,c("lat","long")]
dists <- distances(xy)
q <- samples[,"monticola"]
dat <- data.frame(x=-dists$dist.orig,
                  n=2,
                  q=q)
snps <- read.csv("./GBS/monticola_310_SNP_freqs_rescaled.csv",header=T,row.names=1)
hi <- hi.mle(snps)
loci <- read.csv("./GBS/Locus_Names.csv",header=T,row.names=1)
head(loci)
match(colnames(snps),rownames(loci))
loci$ID[match(colnames(snps),rownames(loci))]
dat$hi <- hi;dat <- data.frame(dat,z=snps)

# Calculating center and width for all clines: same answer as joint_clines
dat$x.s <- trans.down(dat$x)
fit.noTails.all <- function(z){
  dd <- na.omit(data.frame(dat$x.s,dat$n,z))
  colnames(dd) <- c("x.s","n","z")
  return(fit.noTails(dd$x.s,dd$n,dd$z))
}

res.noTails.all <- apply(dat[,5:314],MARGIN=2,FUN=fit.noTails.all)
est.noTails.all <- do.call(rbind.data.frame,res.noTails.all)
est.noTails.all$c <- trans.up(est.noTails.all$c,dat$x)
est.noTails.all$w <- est.noTails.all$w*diff(range(dat$x))

# function for dm testing

dm.test.c <- function(i,j){
  dd.i <- data.frame(na.omit(dat[,c(1:4,i+4,315)]))
  dd.j <- data.frame(na.omit(dat[,c(1:4,j+4,315)]))
  
  # fit clines on 0 to 1
  est.i <- fit.noTails(dd.i$x.s,dd.i$n,dd.i$z)
  est.j <- fit.noTails(dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the swapped center on 0 to 1
  c.swap.j <- trans.down(-trans.up(est.i$c,dat$x),dat$x)
  c.swap.i <- trans.down(-trans.up(est.j$c,dat$x),dat$x)
  
  # calculate the likelihood of no swap
  LL.i <- LL.noTails(unlist(est.i[,1:2]),dd.i$x.s,dd.i$n,dd.i$z)
  LL.j <- LL.noTails(unlist(est.j[,1:2]),dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the likelihood of swap
  LL.dm.i <- LL.noTails(unlist(c(c.swap.i,est.i[,2])),dd.i$x.s,dd.i$n,dd.i$z)
  LL.dm.j <- LL.noTails(unlist(c(c.swap.j,est.j[,2])),dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the likelihood ratio
  LLr.i <- -2*(LL.dm.i-LL.i)
  LLr.j <- -2*(LL.dm.j-LL.j)
  
  # create dataframe
  results <- data.frame(c(i,j),c(est.i[,1],est.j[,1]),
                        c(est.i[,2],est.j[,2]),
                        c(c.swap.i,c.swap.j),c(LL.i,LL.j),
                        c(LL.dm.i,LL.dm.j),c(LLr.i,LLr.j))
  colnames(results) <- c("snps","c","w","dm.c","LL","dm.LL","ratio")
  return(results)
}

# plotting dm.test

dm.plot <- function(i,j){
  
  dm <- dm.test(i,j)
  
  par(mfrow=c(1,1))
  plot(dat$x.s,dat[,4+i],main='Both SNPS',xlab='x.s',ylab='pA',col='chocolate')
  points(dat$x.s,dat[,4+j],col="darkgreen")
  curve(Barton.p.noTails(x,unlist(dm[1,2:3])),col='chocolate',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[2,2:3])),col='darkgreen',add=TRUE)
  
  par(mfrow=c(1,2))
  plot(dat$x.s,dat[,4+i],main=c('snps ', i),xlab='x.s',ylab='pA')
  curve(Barton.p.noTails(x,unlist(dm[1,2:3])),col='blue',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[1,c(4,3)])),col="red",add=TRUE)
  
  plot(dat$x.s,dat[,4+j],main=c('snps ',j),xlab='x.s',ylab='pA')
  curve(Barton.p.noTails(x,unlist(dm[2,2:3])),col='blue',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[2,c(4,3)])),col="red",add=TRUE)
  
}

# testing on width and center for dm testing

dm.test.cw <- function(i,j){
  dd.i <- data.frame(na.omit(dat[,c(1:4,i+4,315)]))
  dd.j <- data.frame(na.omit(dat[,c(1:4,j+4,315)]))
  
  # fit clines on 0 to 1
  est.i <- fit.noTails(dd.i$x.s,dd.i$n,dd.i$z)
  est.j <- fit.noTails(dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the swapped center on 0 to 1
  c.swap.j <- trans.down(-trans.up(est.i$c,dat$x),dat$x)
  c.swap.i <- trans.down(-trans.up(est.j$c,dat$x),dat$x)
  
  # swap the widths
  w.swap.j <- est.i[,2]
  w.swap.i <- est.j[,2]
  
  # calculate the likelihood of no swap
  LL.i <- LL.noTails(unlist(est.i[,1:2]),dd.i$x.s,dd.i$n,dd.i$z)
  LL.j <- LL.noTails(unlist(est.j[,1:2]),dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate likelihood of swapped center only
  LL.dm.i <- LL.noTails(unlist(c(c.swap.i,est.i[,2])),dd.i$x.s,dd.i$n,dd.i$z)
  LL.dm.j <- LL.noTails(unlist(c(c.swap.j,est.j[,2])),dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the likelihood ratio
  LLr.c.i <- -2*(LL.dm.i-LL.i)
  LLr.c.j <- -2*(LL.dm.j-LL.j)
  
  # calculate the likelihood of swapped center and width
  LL.dm.i <- LL.noTails(unlist(c(c.swap.i,w.swap.i)),dd.i$x.s,dd.i$n,dd.i$z)
  LL.dm.j <- LL.noTails(unlist(c(c.swap.j,w.swap.j)),dd.j$x.s,dd.j$n,dd.j$z)
  
  # calculate the likelihood ratio
  LLr.cw.i <- -2*(LL.dm.i-LL.i)
  LLr.cw.j <- -2*(LL.dm.j-LL.j)
  
  
  # create dataframe
  results <- data.frame(c(i,j),c(est.i[,1],est.j[,1]),
                        c(est.i[,2],est.j[,2]),
                        c(c.swap.i,c.swap.j),c(w.swap.i,w.swap.j),
                        c(LLr.c.i,LLr.c.j),c(LLr.cw.i,LLr.cw.j))
  colnames(results) <- c("snps","c","w","swap.c","swap.w","test.c","test.cw")
  return(results)
}

dm.plot.cw <- function(i,j){
  
  dm <- dm.test.cw(i,j)
  
  par(mfrow=c(1,1))
  plot(dat$x.s,dat[,4+i],main='Both SNPS',xlab='x.s',ylab='pA',col='chocolate')
  points(dat$x.s,dat[,4+j],col="darkgreen")
  curve(Barton.p.noTails(x,unlist(dm[1,2:3])),col='chocolate',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[2,2:3])),col='darkgreen',add=TRUE)
  
  par(mfrow=c(1,2))
  plot(dat$x.s,dat[,4+i],main=c('snps ', i),xlab='x.s',ylab='pA')
  curve(Barton.p.noTails(x,unlist(dm[1,2:3])),col='blue',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[1,c(4,3)])),col='orange',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[1,4:5])),col="red",add=TRUE)
  
  legend(.5,.2,legend=c("Barton","C Swap","C&W Swap"),fill=c("blue","orange","red"))
  
  plot(dat$x.s,dat[,4+j],main=c('snps ',j),xlab='x.s',ylab='pA')
  curve(Barton.p.noTails(x,unlist(dm[2,2:3])),col='blue',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[2,c(4,3)])),col='orange',add=TRUE)
  curve(Barton.p.noTails(x,unlist(dm[2,4:5])),col="red",add=TRUE)
}

for (i in 1:5){
  for (j in 6:10){
    if (i<j){
      print(dm.test.cw(i,j))
    }
  }
}

# all the ones coming out as dm have an insignificant width
dm.plot.cw(2,10)
dm.test.cw(2,10)

dm.plot.cw(1,7)
dm.test.cw(1,7)

dm.plot.cw(2,8)
dm.test.cw(2,8)

# filtering for significant width
est.noTails.all$snps <- seq(1,310,1)
est.noTails.all$Pw <- 1-pnorm(est.noTails.all$w,est.admixture$w,est.admixture$w/3)
est.sig.w <- filter(est.noTails.all,(Pw>.05) & (c>0))
est.sig.w

# finding all with significant width

k=1
dm.snps <- data.frame(rep(NA,10000),rep(NA,10000))
for (i in 1:310){
  for (j in 1:310){
    if (i<j) {
      dm <- dm.test.cw(i,j)
      if (all(dm[,6] < 3.841) & all(dm[,7] < 5.991)){
        dm.snps[k,1] <- i
        dm.snps[k,2] <- j
        print(c(i,j))
        k=k+1
      }
    }
  }
}
dm.snps <- na.omit(dm.snps)
colnames(dm.snps) <- c("snps1","snps2")

dm.test.all <- data.frame(dm.snps)

for (i in 1:nrow(dm.snps)){
  print(dm.test.cw(dm.test.all[i,1],dm.test.all[i,2]))
}

################################################################################
####################### Fitting environmental cline ############################
################################################################################

# FUNCTIONS 
############

# joint environmental and genomic functions
pars.J.e <- data.frame(c=runif(1000),
                       w=runif(1000,0,2),
                       mu=runif(1000,-100,100),
                       nu=runif(1000,0,1000))

LL.Joint.e.noTails <- function(pars.J.e,x.e,x.g,n,y)
{p <- Joint.p.noTails(x.e,x.g,pars.J.e)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

fitJoint.e.noTails <- function(pars.J.e,x.e,x.g,n,y)
{res.Joint <- multistart(par=pars.J.e,fn=LL.Joint.e.noTails,
                         x.e=x.e,x.g=x.g,n=n,y=y,
                         lower=c(-0.1,0,-Inf,0),
                         upper=c(1.1,2.2,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint[which.max(res.Joint$value),],2)
est.J}

# joint environmental, spatial, and genomic cline functions

pars.J.3 <- data.frame(c.e=runif(5000),
                       w.e=runif(5000,0,2),
                       c.x=runif(5000),
                       c.x=runif(5000,0,2),
                       mu=runif(5000,-100,100),
                       nu=runif(5000,0,1000))

Joint.3p <- function(x.e,x.s,x.g,pars.J.3){
  Barton.p.noTails(x.e,pars.J.3[1:2])*Barton.p.noTails(x.s,pars.J.3[3:4])*Fitz.p(x.g,pars.J.3[5:6])}

LL.Joint.3 <- function(pars.J.3,x.e,x.s,x.g,n,y)
{p <- Joint.3p(x.e,x.s,x.g,pars.J.3)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

fitJoint.3 <- function(pars.J.3,x.e,x.s,x.g,n,y)
{res.Joint <- multistart(par=pars.J.3,fn=LL.Joint.3,
                         x.e=x.e,x.s=x.s,x.g=x.g,n=n,y=y,
                         lower=c(-0.1,0,-.1,0,-Inf,0),
                         upper=c(1.1,2.2,1.1,2.2,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint[which.max(res.Joint$value),],2)
est.J}

# multistart likelihoods
res.1 <- function(x.s,n,y)
{res.noTails <- multistart(parmat=pars.N,
                           fn=LL.noTails,
                           x=x.s,n=n,y=y,
                           lower=c(-0.1,0),
                           upper=c(1.1,2.2),
                           method="L-BFGS-B",
                           control=list(maximize=TRUE))
est.noTails <- round(res.noTails,2)
est.noTails
}

res.2 <- function(pars.J,x.s,x.g,n,y)
{res.Joint <- multistart(par=pars.J,fn=LL.Joint.noTails,
                         x.s=x.s,x.g=x.g,n=n,y=y,
                         lower=c(-0.1,0,-Inf,0),
                         upper=c(1.1,2.2,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint,2)
est.J}

res.3 <- function(pars.J.3,x.e,x.s,x.g,n,y)
{res.Joint <- multistart(par=pars.J.3,fn=LL.Joint.3,
                         x.e=x.e,x.s=x.s,x.g=x.g,n=n,y=y,
                         lower=c(-0.1,0,-.1,0,-Inf,0),
                         upper=c(1.1,2.2,1.1,2.2,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint,2)
est.J}


# DATA 
##########

#Functions
source("./aux_funcs.R")

#sample data for monticola - 71 samples from Pyron et al. 2023
samples <- read.csv("./GBS/gbs_71_localities.csv",header=T,row.names=1)

#extract lat/longs
xy <- samples[,c("lat","long")]

#get distances from furthest point
#you may wish to provide your own.
#This returns a matrix and a vector
dists <- distances(xy)

#extract individual ancestry coefficients
q <- samples[,"monticola"]

# extract the environmental covariate
e <- samples[,"BIO15"]

#make initial data object
dat <- data.frame(x=-dists$dist.orig,
                  n=2,
                  q=q,
                  e=e)
# normalize x
dat$x.s <- trans.down(dat$x)

#SNPs for monticola - 310 loci from Pyron et al. 2023, rescaled to monticola
snps <- read.csv("./GBS/monticola_310_SNP_freqs_rescaled.csv",header=T,row.names=1)

#Hybrid index for sampled loci, assuming biallelic SNPs
hi <- hi.mle(snps)

# normalize and reverse environmental 
dat$x.e <- (dat$e - min(dat$e)) / diff(range(dat$e))
dat$x.e <- 1 - dat$x.e                            # not applicable in all cases

#Locus names - to be used later for locus-specific estimates
loci <- read.csv("./GBS/Locus_Names.csv",header=T,row.names=1)
head(loci)

#Match the loci to our SNPs
match(colnames(snps),rownames(loci))
loci$ID[match(colnames(snps),rownames(loci))]

#Add HI and SNPs to data object
dat$hi <- hi;dat <- data.frame(dat,z=snps)


# FITTING THE CLINES
####################

# top 10 likelihoods for snps 1-50

ridge.test <- list()

for (i in 1:50){
  
  dd <- data.frame(na.omit(dat[,c(1:7,i+7)]))
  
  res <- res.3(pars.J.3,dd$x.e,dd$x.s,dd$q,dd$n,dd$z)
  res <- res[order(res$value,decreasing=TRUE),]
  
  ridge.test <- append(ridge.test, list(res[1:10,]))
  
}

ridge.test
ridges

# plotting the top 10 snps
j <- 15
dd <- data.frame(na.omit(dat[,c(1:7,j+7)]))
res <- data.frame(ridge.test[j])

for (i in 1:10){
  
  par(mfrow=c(3,1))
  
  plot(dd$x.e,dd$z,main=c("Joint: Environmental",i))
  curve(Barton.p.noTails(x,unlist(res[i,1:2])),from=0,to=1,add=TRUE) 
  plot(dd$x.s,dd$z,main="Joint: Spatial")
  curve(Barton.p.noTails(x,unlist(res[i,3:4])),from=0,to=1,add=TRUE)
  plot(dd$q,dd$z,main="Joint: Genomic")
  curve(Fitz.p(x,unlist(res[i,5:6])),from=0,to=1,add=TRUE)
  
}

## FOR INVESTIGATION PURPOSES
# environmental cline only
est.e <- fit.noTails(dd$x.e,dd$n,dd$z)

par(mfrow=c(1,1))
plot(dd$x.e,dd$z,main="Environmental Covariate",xlab="BIO15",ylab="pA")
curve(Barton.p.noTails(x,unlist(est.e[,1:2])),from=0,to=1,add=TRUE)

res <- res.1(dd$x.e,dd$n,dd$z) # this worked, optimized to the same parameters every time
res <- res[order(res$value),]

# estimating the joint parameters for spatial and genomic
est.j <- fitJoint.noTails(pars.J,dd$x.s,dd$q,dd$n,dd$z)

par(mfrow=c(2,1))
plot(dd$x,dd$z,main="Joint Spatial")
curve(Barton.p.noTails(x,unlist(est.j[,1:2])),from=min(dd$x),to=max(dd$x),add=TRUE)
plot(dd$q,dd$z,main="Joint Genomic")
curve(Fitz.p(x,unlist(est.j[,3:4])),from=0,to=1,add=TRUE)

res <- res.2(pars.J,dd$x.s,dd$q,dd$n,dd$z) # this worked, parameters all close and differing likelihoods
res <- res[order(res$value,decreasing=TRUE),]

## estimating joint parameters for environmental and genomic
est.e <- fitJoint.e.noTails(pars.J.e,dd$x.e,dd$q,dd$n,dd$z)

plot(dd$x.e,dd$z,main="Joint Environmental Covariate")
curve(Barton.p.noTails(x,unlist(est.e[,1:2])),from=0,to=1,add=TRUE) 
plot(dd$q,dd$z,main="Joint Genomic")
curve(Fitz.p(x,unlist(est.e[,3:4])),from=0,to=1,add=TRUE)

res <- res.2(pars.J.e,dd$x.e,dd$q,dd$n,dd$z) 
res <- res[order(res$value,decreasing=TRUE),]

# MODEL SELECTION / VARIABLE IMPORTANCE 
#######################################

# making the joint spatial and environmental

pars.J.es <- data.frame(c.e=runif(1000),
                        w.e=runif(1000,0,2),
                        c.s=runif(1000),
                        w.s=runif(1000,0,2))

Joint.p.es <- function(x.e,x.s,pars.J.es){Barton.p.noTails(x.e,pars.J.es[1:2])*Barton.p(x.s,pars.J.es[3:4])}

LL.Joint.es <- function(pars.J.es,x.e,x.s,n,y)
{p <- Joint.p.noTails(x.e,x.s,pars.J.es)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

fitJoint.es <- function(pars.J.es,x.e,x.s,n,y)
{res.Joint <- multistart(par=pars.J.es,fn=LL.Joint.es,
                         x.s=x.s,x.e=x.e,n=n,y=y,
                         lower=c(-0.1,0,-.1,0),
                         upper=c(1.1,2.2,1.1,2.2),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint[which.max(res.Joint$value),],2)
est.J}

# loop through each model, calculate AIC, rank by AIC

best.model <- function(x.s,x.e,x.g,n,z){
  
  est1 <- fit.noTails(x.s,n,z) # c.s, w.s
  colnames(est1) <- c("c.s","w.s","value","fevals","gevals","convergence")
  
  est2 <- fit.noTails(x.e,n,z) # c.e, w.e
  colnames(est2) <- c("c.e","w.e","value","fevals","gevals","convergence")
  
  est3 <- fitFitz(x.g,n,z)
  
  est4 <- fitJoint.noTails(pars.J,x.s,x.g,n,z) # c.s, w.s, mu, nu
  colnames(est4) <- c("c.s","w.s","mu","nu","value","fevals","gevals","convergence")
  
  est5 <- fitJoint.e.noTails(pars.J.e,x.e,x.g,n,z) # c.e, w.e, mu, nu
  colnames(est5) <- c("c.e","w.e","mu","nu","value","fevals","gevals","convergence")
  
  est6 <-fitJoint.es(pars.J.es,x.e,x.s,n,z)
  
  est7 <- fitJoint.3(pars.J.3,x.e,x.s,x.g,n,z) # c.e, w.e, c.s, w.s, mu, nu
  colnames(est7) <- c("c.e","w.e","c.s","w.s","mu","nu","value","fevals","gevals","convergence")
  
  results <- rbind.fill(est1,est2,est3,est4,est5,est6,est7)
  results <- results[,c(1,2,7,8,9,10,3,4,5,6)]
  results$k <- c(1,1,1,2,2,2,3)
  results$AIC <- (2*results$k) - (2*results$value)
  return(results[order(results$AIC),])
}

dd <- data.frame(na.omit(dat[,c(1:7,2+7)]))
results <- best.model(dd$x.s,dd$x.e,dd$q,dd$n,dd$z)

# modeling results (need to change the order for each)
par(mfrow=c(3,1))
plot(dd$x.s,dd$z,main="Spatial")
curve(Barton.p.noTails(x,unlist(results[7,1:2])),from=0,to=1,add=TRUE)
plot(dd$x.e,dd$z,main="Environmental")
curve(Barton.p.noTails(x,unlist(results[5,3:4])),from=0,to=1,add=TRUE)
plot(dd$q,dd$z,main="Genomic")
curve(Fitz.p(x,unlist(results[6,5:6])),from=0,to=1,add=TRUE)

par(mfrow=c(2,1))
plot(dd$x.s,dd$z,main="Joint Spatial")
curve(Barton.p.noTails(x,unlist(results[4,1:2])),from=0,to=1,add=TRUE)
plot(dd$x.e,dd$z,main="Joint Environmental")
curve(Barton.p.noTails(x,unlist(results[4,3:4])),from=0,to=1,add=TRUE)

par(mfrow=c(2,1))
plot(dd$x.s,dd$z,main="Joint Spatial")
curve(Barton.p.noTails(x,unlist(results[2,1:2])),from=0,to=1,add=TRUE)
plot(dd$q,dd$z,main="Joint Genomic")
curve(Fitz.p(x,unlist(results[2,5:6])),from=0,to=1,add=TRUE)

par(mfrow=c(2,1))
plot(dd$x.e,dd$z,main="Joint Environmental")
curve(Barton.p.noTails(x,unlist(results[1,3:4])),from=0,to=1,add=TRUE)
plot(dd$q,dd$z,main="Joint Genomic")
curve(Fitz.p(x,unlist(results[1,5:6])),from=0,to=1,add=TRUE)

par(mfrow=c(3,1))
plot(dd$x.s,dd$z,main="Joint Spatial")
curve(Barton.p.noTails(x,unlist(results[3,1:2])),from=0,to=1,add=TRUE)
plot(dd$x.e,dd$z,main="Joint Environmental")
curve(Barton.p.noTails(x,unlist(results[3,3:4])),from=0,to=1,add=TRUE)
plot(dd$q,dd$z,main="Joint Genomic")
curve(Fitz.p(x,unlist(results[3,5:6])),from=0,to=1,add=TRUE)

################################################################################
############### REBUILDING ENVIRONMENTAL CLINE IN LINEAR SPACE #################
################################################################################

# DATA
######

#Functions
source("./aux_funcs.R")

#sample data for monticola - 71 samples from Pyron et al. 2023
samples <- read.csv("./GBS/gbs_71_localities.csv",header=T,row.names=1)

#extract lat/longs
xy <- samples[,c("lat","long")]

#get distances from furthest point
#you may wish to provide your own.
#This returns a matrix and a vector
dists <- distances(xy)

#extract individual ancestry coefficients
q <- samples[,"monticola"]

# extract the environmental covariate
e <- samples[,"BIO15"]

#make initial data object
dat <- data.frame(x=-dists$dist.orig,
                  n=2,
                  q=q,
                  e=e)
# normalize x
dat$x.s <- trans.down(dat$x)

#SNPs for monticola - 310 loci from Pyron et al. 2023, rescaled to monticola
snps <- read.csv("./GBS/monticola_310_SNP_freqs_rescaled.csv",header=T,row.names=1)

#Hybrid index for sampled loci, assuming biallelic SNPs
hi <- hi.mle(snps)

# normalize and reverse environmental 
dat$x.e <- (dat$e - min(dat$e)) / diff(range(dat$e))
dat$x.e <- 1 - dat$x.e                            # not applicable in all cases

#Locus names - to be used later for locus-specific estimates
loci <- read.csv("./GBS/Locus_Names.csv",header=T,row.names=1)
head(loci)

#Match the loci to our SNPs
match(colnames(snps),rownames(loci))
loci$ID[match(colnames(snps),rownames(loci))]

#Add HI and SNPs to data object
dat$hi <- hi;dat <- data.frame(dat,z=snps)


# updating some values to make my life easier
dd <- data.frame(na.omit(dat[,c(1:7,3+7)]))
dd$z[dd$z==0] <- inv.logit(-5)
dd$z[dd$z==1] <- inv.logit(5)
dd$q[dd$q==0] <- inv.logit(-5)
dd$q[dd$q==1] <- inv.logit(5)
dd$x.s[dd$x.s==0] <- inv.logit(-5)
dd$x.s[dd$x.s==1] <- inv.logit(5)


# logit function
logit <- function(x){log(x/(1-x))}

# Barton spatial logit function
spatial.logit <- function(x.s,pars){pars[2]*(x.s-pars[1])}

# genomic logit function
genomic.logit <- function(logit.x.g,pars){pars[2]*logit.x.g - pars[1]}


# our normal spatial cline
est <- fit.noTails(dd$x.s,dd$n,dd$z)

# plotting spatial cline in logit space
plot(dd$x.s,logit(dd$z),main="fitting estimated parameters to logit")
curve(spatial.logit(x,unlist(c(est[,1],4/est[,2]))),from=min(dd$x.s),to=max(dd$x.s),add=TRUE,col="blue")

# just calculating p the normal barton way and then taking the logit
p <- Barton.p.noTails(dd$x.s,unlist(est[,1:2]))
p[p==0] <- .01
p[p==1] <- .99
p.logit <- logit(p)
points(dd$x.s,p.logit, col="red")

# trying this for genomic: apparently the line in genomic is just equal to nu*logit(q) - mu
est.g <- fitFitz(dd$q,dd$n,dd$z)
p <- Fitz.p(dd$q,unlist(est.g[,1:2]))
logit.p <- logit(p)
plot(logit(dd$q),logit(dd$z))
points(logit(dd$q),logit.p,col="red")
curve(genomic.logit(x,unlist(est.g[,1:2])),from=min(logit(dd$q)),to=max(logit(dd$q)),add=TRUE)

# lets try to fit these functions now, starting with the spatial
pars.S <- data.frame(c=runif(100),
                     w=runif(100,.0001,2))

fitlogit.s <- function(x.s,n,y)
{res.logit <- multistart(parmat=pars.S,
                         fn=LL.noTails,
                         x=x.s,n=n,y=y,
                         lower=c(-0.1,0.0001),
                         upper=c(1.1,2.2),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.logit <- round(res.logit[which.max(res.logit$value),],2)
est.logit$b <- 4/est.logit$w
est.logit <- est.logit[,c(1,7,3:5)]
est.logit}

est.s <- fitlogit.s(dd$x.s,dd$n,dd$z)
plot(dd$x.s,logit(dd$z))
curve(spatial.logit(x,unlist(est.s[,1:2])),from=0,to=1,add=TRUE)

# plotting genomic in logit space
est.g <- fitFitz(dd$q,dd$n,dd$z)
p <- Fitz.p(dd$q,unlist(est.g[,1:2]))
logit.p <- logit(p)
plot(logit(dd$q),logit(dd$z))
lines(logit(dd$q),logit.p,col="red")
curve(genomic.logit(x,unlist(est.g[,1:2])),from=min(logit(dd$q)),to=max(logit(dd$q)),add=TRUE)

################################################################################
################# JOINT CLINE ESTIMATION IN LOGIT SPACE ########################
################################################################################

pars.SG <- pars.J <- data.frame(c=runif(1000),
                                w=runif(1000,0.0001,2),
                                mu=runif(1000,-100,100),
                                nu=runif(1000,0,1000))

joint.logit <- function(x.s,logit.x.g,pars.SG){((4/pars.SG[2])*(x.s-pars.SG[1])) + (pars.SG[4]*logit.x.g - pars.SG[3])}

fitlogit.sg <- function(x.s,x.g,n,y)
{res.Joint <- multistart(par=pars.J,fn=LL.Joint.noTails,
                         x.s=x.s,x.g=x.g,n=n,y=y,
                         lower=c(-0.1,0,-Inf,0),
                         upper=c(1.1,2.2,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint[which.max(res.Joint$value),],2)
est.J$b <- 4/est.J$w
est.J <- est.J[,c(1,9,3:8)]}

# plotting the joint cline in logit space. I added the logit(p) to verify these are correct.
dd <- data.frame(na.omit(dat[,c(1:7,14+7)]))
dd$z[dd$z==0] <- inv.logit(-5)
dd$z[dd$z==1] <- inv.logit(5)
dd$q[dd$q==0] <- inv.logit(-5)
dd$q[dd$q==1] <- inv.logit(5)
dd$x.s[dd$x.s==0] <- inv.logit(-5)
dd$x.s[dd$x.s==1] <- inv.logit(5)

est.J <- fitJoint.noTails(pars.J,dd$x.s,dd$q,dd$n,dd$z) # normal parameters
est.J$c <- trans.down(est.J$c,dat$x)
est.J$w <- est.J$w / diff(range(dat$x))
est.J.logit <- fitlogit.sg(dd$x.s,dd$q,dd$n,dd$z) # logit space parameters: these check out the same

par(mfrow=c(2,2))
plot(dd$x.s,logit(dd$z),main="Joint spatial in logit space")
curve(spatial.logit(x,unlist(est.J.logit[,1:2])),from=min(dd$x.s),to=max(dd$x.s),add=TRUE)
p <- Barton.p.noTails(dd$x.s,unlist(est.J[,1:2]))
lines(dd$x.s,logit(p),col="deeppink1")
abline(spatial_model,col="blue")

plot(dd$x.s,dd$z,main="Joint spatial")
curve(Barton.p.noTails(x,unlist(est.J[,1:2])),from=0,to=1,add=TRUE)
curve(Barton.p.noTails(x,unlist(spatial[,1:2])),from=0,to=1,add=TRUE,col="blue")

p <- Fitz.p(dd$q,unlist(est.J[,3:4]))
logit.p <- logit(p)
plot(logit(dd$q),logit(dd$z),main="Genomic spatial in logit space")
curve(genomic.logit(x,unlist(est.J.logit[,3:4])),from=min(logit(dd$q)),to=max(logit(dd$q)),add=TRUE)
lines(logit(dd$q),logit.p,col="deeppink1")
abline(-mu,nu,col="blue")

plot(dd$q,dd$z,main="Joint genomic")
curve(Fitz.p(x,unlist(est.J[,3:4])),from=0,to=1,add=TRUE)
curve(Fitz.p(x,unlist(genomic[,1:2])),from=0,to=1,add=TRUE,col="blue")

# LOGIT OPTIMIZATION
####################

logit <- function(x){log(x/(1-x))}
inv.logit <- function(x){exp(x)/(1 + exp(x))}

pars.L <- data.frame(c = runif(100),
                     b = runif(100,2,400))

Barton.p.logit <- function(x.s,pars.L){pars.L[2]*(x.s-pars.L[1])}

LL.noTails.logit <- function(x.s,y,n,pars.L)
{ p <- inv.logit(Barton.p.logit(x.s,pars.L))
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

fit.noTails.logit <- function(x.s,n,y,pars.L)
{res.noTails <- multistart(parmat=pars.L,
                           fn=LL.noTails.logit,
                           x.s=x.s,n=n,y=y,
                           lower=c(-0.1,2),
                           upper=c(1.1,400),
                           method="L-BFGS-B",
                           control=list(maximize=TRUE))
est.noTails <- round(res.noTails[which.max(res.noTails$value),],2)
est.noTails}

est.logit <- fit.noTails.logit(dd$x.s,dd$n,dd$z,pars.L)
est <- fit.noTails(dd$x.s,dd$n,dd$z)

# check equivalence
par(mfrow=c(1,1))
plot(dd$x.s,logit(dd$z))
curve(Barton.p.logit(x,unlist(c(est[,1],4/est[,2]))),from=0,to=1,add=TRUE,col="red")  
curve(Barton.p.logit(x,unlist(est.logit[,1:2])),from=0,to=1,add=TRUE,col="blue")  

# Optimizing for the genomic cline in logit space
Fitz.logit.p <- function(logit.x.g,pars.G){pars.G[2]*logit.x.g - pars.G[1]}

LL.Fitz.logit <- function(logit.x.g,y,n,pars.G){
  p <- inv.logit(Fitz.logit.p(logit.x.g,pars.G))
  p <- replace(p, p <= 0, .Machine$double.xmin)
  p <- replace(p, p >= 1, 1 - .Machine$double.neg.eps)
  LL <- likelihood.C(p,y,n)
  return(sum(LL))
}

fit.Fitz.logit <- function(logit.x.g,y,n,pars.G){
  res.Fitz <- multistart(par=pars.G,
                         fn=LL.Fitz.logit,
                         logit.x.g=logit.x.g,n=n,y=y,
                         lower=c(-Inf,0.0),
                         upper=c(Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
  est.Fitz <- round(res.Fitz[which.max(res.Fitz$value),],2)
  est.Fitz
}

# Checking for equivalence
est.logit <- fit.Fitz.logit(logit(dd$q),dd$z,dd$n,pars.G)
est <- fitFitz(dd$q,dd$n,dd$z)

plot(logit(dd$q),logit(dd$z),main="Genomic fit in logit space")
abline(-est[,1],est[,2],col="red")
abline(-est.logit[,1],est.logit[,2],col="blue")

# Optimizing the joint cline model in logistic space
pars.JL <- data.frame(c = runif(100),
                      b = runif(100,2,400),
                      mu=runif(100,-100,100),
                      nu=runif(100,0,1000))

Joint.p.logit <- function(x.s,logit.x.g,pars.JL){pars.JL[2]*(x.s-pars.JL[1]) + (pars.JL[3]*logit.x.g - pars.JL[4])}

LL.Joint.logit <- function(pars.JL,x.s,logit.x.g,n,y)
{p <- inv.logit(Joint.p.logit(x.s,logit.x.g,pars.JL))
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

fitJoint.logit <- function(pars.JL,x.s,logit.x.g,n,y)
{res.Joint <- multistart(par=pars.JL,fn=LL.Joint.logit,
                         x.s=x.s,logit.x.g=logit.x.g,n=n,y=y,
                         lower=c(-0.1,2,-Inf,0),
                         upper=c(1.1,400,Inf,Inf),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
est.J <- round(res.Joint[which.max(res.Joint$value),],2)
est.J}


fitJoint.logit(pars.JL,dd$x.s,logit(dd$q),dd$n,dd$z)
fitJoint.noTails(pars.J,dd$x.s,dd$q,dd$n,dd$z)

plot(dd$x.s,dd$z)

# Visualizing Barton no tails approach
################################################################################
source("./aux_funcs.R")
source("./marlee_aux.R")

est.noTails <- fit.noTails(dd$x.s,dd$n,dd$z)
est.noTails     

est.noTails.minmax <- fit.noTails.minmax(dd$x.s,dd$n,dd$z)
est.noTails.minmax

plot(dd$x.s,dd$z/2)
curve(Barton.p.noTails.minmax(x,unlist(est.noTails.minmax[,1:4])),from=0,to=1,add=TRUE)

for (i in 1:310){
  dd <- data.frame(na.omit(dd.all[,c(1:5,i+5)]))
  
  est.noTails <- fit.noTails(dd$x.s,dd$n,dd$z)
  est.noTails.minmax <- fit.noTails.minmax(dd$x.s,dd$n,dd$z)
  
  plot(dd$x.s,dd$z/2,xlab="x.s",ylab="pA",main=i,col="black",pch=20,cex=2)
  curve(Barton.p.noTails(x,unlist(est.noTails[,1:2])),from=0,to=1,add=TRUE,col="blue",lwd=1.5)
  curve(Barton.p.noTails.minmax(x,unlist(est.noTails.minmax[,1:4])),from=0,to=1,add=TRUE,col="red",lwd=1.5)
  readline(prompt="Press [enter] to continue...")
}

# TAILS OF JOINT SPATIAL/GENOMIC LOGISTIC REGRESSION FIT
################################################################################

# admixture cline
est.joint <- fit.bTails.Joint(dd$x.s,dd$x.g,dd$n,dd$y)
est.Barton <- fit.bTails(dd$x.s,dd$n,dd$y)
est.Fitz <- fitFitz(dd$x.g,dd$n,dd$y)

plot(dd$x.s,dd$y,xlab="x.s",ylab="pA",col="black",pch=20,cex=2)
curve(Barton.p.bTails(x,unlist(est.joint[,1:6])),from=0,to=1,add=TRUE,col="red")
curve(Barton.p.bTails(x,unlist(est.Barton[,1:6])),from=0,to=1,add=TRUE,col="blue")

plot(dd$x.g,dd$y,xlab="x.g",ylab="pA",col="black",pch=20,cex=2)
curve(Fitz.p(x,unlist(est.joint[,7:8])),from=0,to=1,add=TRUE,col="red")
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,add=TRUE,col="blue")

# Inidividual snps
for (i in 1:310){
  dd <- data.frame(na.omit(dd.all[,c(1:5,5+i)]))
  
  est.joint <- fit.bTails.Joint(dd$x.s,dd$x.g,dd$n,dd$z/2)
  est.Barton <- fit.bTails(dd$x.s,dd$n,dd$z)
  est.Fitz <- fitFitz(dd$x.g,dd$n,dd$z)
  
  par(mfrow=c(2,1))
  plot(dd$x.s,dd$z/2,xlab="x.s",ylab="pA",col="black",pch=20,cex=2,main="Spatial")
  curve(Barton.p.bTails(x,unlist(est.joint[,1:6])),from=0,to=1,add=TRUE,col="red")
  curve(Barton.p.bTails(x,unlist(est.Barton[,1:6])),from=0,to=1,add=TRUE,col="blue")
  
  plot(dd$x.g,dd$z/2,xlab="x.g",ylab="pA",col="black",pch=20,cex=2,main="Genomic")
  curve(Fitz.p(x,unlist(est.joint[,7:8])),from=0,to=1,add=TRUE,col="red")
  curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,add=TRUE,col="blue")
  mtext(text=i,outer=TRUE,line=-1)
  
  readline(prompt="Press [enter] to continue...")
}

# Trying each of the logistic tail options
# No tails, both tails, left tail, right tail, mixed
i=11
dd <- data.frame(na.omit(dd.all[,c(1:5,5+i)]))

est.No <- fitJoint.noTails(pars.J,dd$x.s,dd$x.g,dd$n,dd$z/2)
est.B <- fit.bTails.Joint(dd$x.s,dd$x.g,dd$n,dd$z/2)
est.M <- fit.mTails.Joint(dd$x.s,dd$x.g,dd$n,dd$z/2)
est.L <- fit.lTails.Joint(dd$x.s,dd$x.g,dd$n,dd$z/2)
est.R <- fit.rTails.Joint(dd$x.s,dd$x.g,dd$n,dd$z/2)

est.Barton <- fit.bTails(dd$x.s,dd$n,dd$z)
est.Fitz <- fitFitz(dd$x.g,dd$n,dd$z)

par(mfrow=c(2,1))
plot(dd$x.s,dd$z/2,main="Spatial",pch=20)
curve(Barton.p.noTails(x,unlist(est.Barton[,1:2])),from=0,to=1,add=TRUE,col="black")
curve(Barton.p.noTails(x,unlist(est.No[,1:2])),from=0,to=1,add=TRUE,col="purple")
curve(Barton.p.bTails(x,unlist(est.B[,1:6])),from=0,to=1,add=TRUE,col="blue")
curve(Barton.p.mTails(x,unlist(est.M[,1:4])),from=0,to=1,add=TRUE,col="orange")
curve(Barton.p.lTails(x,unlist(est.L[,1:4])),from=0,to=1,add=TRUE,col="red")
curve(Barton.p.rTails(x,unlist(est.R[,1:4])),from=0,to=1,add=TRUE,col="green")
legend("bottomright",legend=c("Barton","J No Tails","J Both","J Mirror","J Left","J Right"),
       col=c("black","purple","blue","orange","red","green"),lty=1)

plot(dd$x.g,dd$z/2,main="Genomic",pch=20)
curve(Fitz.p(x,unlist(est.Fitz[,1:2])),from=0,to=1,add=TRUE,col="black")
curve(Fitz.p(x,unlist(est.No[,3:4])),from=0,to=1,add=TRUE,col="purple")
curve(Fitz.p(x,unlist(est.B[,7:8])),from=0,to=1,add=TRUE,col="blue")
curve(Fitz.p(x,unlist(est.M[,5:6])),from=0,to=1,add=TRUE,col="orange")
curve(Fitz.p(x,unlist(est.L[,5:6])),from=0,to=1,add=TRUE,col="red")
curve(Fitz.p(x,unlist(est.R[,5:6])),from=0,to=1,add=TRUE,col="green")
legend("bottomright",legend=c("Fitz","J No Tails","J Both","J Mirror","J Left","J Right"),
       col=c("black","purple","blue","orange","red","green"),lty=1)

# Admixture cline
est.No.adm <- fitJoint.noTails(pars.J,dd$x.s,dd$x.g,dd$n,dd$y)
est.B.adm <- fit.bTails.Joint(dd$x.s,dd$x.g,dd$n,dd$y)
est.M.adm <- fit.mTails.Joint(dd$x.s,dd$x.g,dd$n,dd$y)
est.L.adm <- fit.lTails.Joint(dd$x.s,dd$x.g,dd$n,dd$y)
est.R.adm <- fit.rTails.Joint(dd$x.s,dd$x.g,dd$n,dd$y)

est.Barton.adm <- fit.bTails(dd$x.s,dd$n,dd$y)
est.Fitz.adm <- fitFitz(dd$x.g,dd$n,dd$y)

par(mfrow=c(2,1))
plot(dd$x.s,dd$y,main="Spatial")
curve(Barton.p.noTails(x,unlist(est.Barton.adm[,1:2])),from=0,to=1,add=TRUE,col="black")
curve(Barton.p.noTails(x,unlist(est.No.adm[,1:2])),from=0,to=1,add=TRUE,col="purple")
curve(Barton.p.bTails(x,unlist(est.B.adm[,1:6])),from=0,to=1,add=TRUE,col="blue")
curve(Barton.p.mTails(x,unlist(est.M.adm[,1:4])),from=0,to=1,add=TRUE,col="orange")
curve(Barton.p.lTails(x,unlist(est.L.adm[,1:4])),from=0,to=1,add=TRUE,col="red")
curve(Barton.p.rTails(x,unlist(est.R.adm[,1:4])),from=0,to=1,add=TRUE,col="green")
legend("bottomright",legend=c("Barton","J No Tails","J Both","J Mirror","J Left","J Right"),
       col=c("black","blue","orange","red","green"),lty=1)

plot(dd$x.g,dd$y,main="Genomic")
curve(Fitz.p(x,unlist(est.Fitz.adm[,1:2])),from=0,to=1,add=TRUE,col="black")
curve(Fitz.p(x,unlist(est.B.adm[,7:8])),from=0,to=1,add=TRUE,col="blue")
curve(Fitz.p(x,unlist(est.M.adm[,5:6])),from=0,to=1,add=TRUE,col="orange")
curve(Fitz.p(x,unlist(est.L.adm[,5:6])),from=0,to=1,add=TRUE,col="red")
curve(Fitz.p(x,unlist(est.R.adm[,5:6])),from=0,to=1,add=TRUE,col="green")

# CALCULATING AIC FOR ALL TAIL VARIATIONS AND COMPARING
#######################################################
dd <- data.frame(na.omit(dd.all[,c(1:5,5+1)]))
est1 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+2)]))
est2 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+3)]))
est3 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+4)]))
est4 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+5)]))
est5 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+8)]))
est8 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+11)]))
est11 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+43)]))
est43 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

dd <- data.frame(na.omit(dd.all[,c(1:5,5+45)]))
est45 <- fitTails.Joint(dd$x,dd$x.g,dd$n,dd$z)

