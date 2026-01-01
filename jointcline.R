#Barton spatial admixture cline, untransformed
pars.A.minmax  <- data.frame(c=runif(100,0,1),
                             w=runif(100,0,2),
                             p0=runif(100,0,1),
                             p1=runif(100,0,1))

pars.JB <- data.frame(c=runif(500,0,1),
                      w=runif(500,0,2),
                      dL=runif(500,0,2),
                      tL=runif(500),
                      dR=runif(500,0,2),
                      tR=runif(500),
                      mu=runif(500,-100,100),
                      nu=runif(500,0,1000))

pars.JM <- data.frame(c=runif(500),
                     w=runif(500,0,2),
                     d=runif(500,0,2),
                     t=runif(500),
                     mu=runif(500,-100,100),
                     nu=runif(500,0,1000))

pars.JL <- data.frame(c=runif(500),
           w=runif(500,0,2),
           d=runif(500,0,2),
           t=runif(500),
           mu=runif(500,-100,100),
           nu=runif(500,0,1000))

pars.JR <- data.frame(c=runif(500),
                      w=runif(500,0,2),
                      d=runif(500,0,2),
                      t=runif(500),
                      mu=runif(500,-100,100),
                      nu=runif(500,0,1000))

### Probability densities
#PDF for the Barton spatial cline, no tails
Barton.p.noTails.minmax <- function(x,pars){pars[3] + (pars[4] - pars[3])*(1/(1 + exp(-((x - pars[1]) * 4/pars[2]))))}

# PDF for joint cline, both tails
Joint.p.bTails <- function(x,x.g,pars.JB){Barton.p.bTails(x,pars.JB[1:6])*Fitz.p(x.g,pars.JB[7:8])}

# PDF for joint cline, mirror tails
Joint.p.mTails <- function(x,x.g,pars.JM){Barton.p.mTails(x,pars.JM[1:4])*Fitz.p(x.g,pars.JM[5:6])}

#PDF for the joint cline, left tail
Joint.p.lTails <- function(x,x.g,pars.JL){Barton.p.lTails(x,pars.JL[1:4])*Fitz.p(x.g,pars.JL[5:6])}

#PDF for the joint cline, right tail
Joint.p.rTails <- function(x,x.g,pars.JR){Barton.p.rTails(x,pars.JR[1:4])*Fitz.p(x.g,pars.JR[5:6])}


### Likelihoods

#Likelihood for continuous frequencies
likelihood.C <- function(p,y,n)
{#Piecewise
  c1 <- y==0
  c2 <- y==1
  c3 <- y>0 & y<1
  
  #Likelihoods
  lnL <- numeric(length(y))
  lnL[c1] <- log(1-p[c1])*n[c1]
  lnL[c2] <- log(p[c2])*n[c2]
  lnL[c3] <- ((1-y[c3])*log((1-p[c3])/(1-y[c3]))+y[c3]*log(p[c3]/y[c3]))*n[c3]
  
  #Results
  return(lnL)}

#Likelihood for spatial cline, no tails
LL.noTails.minmax <- function(pars.A.minmax,x,n,y)
{p <- Barton.p.noTails.minmax(x,pars.A.minmax)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for joint cline, both tails
LL.bTails.Joint <- function(pars.JB,x,x.g,n,y)
{p <- Joint.p.bTails(x,x.g,pars.JB)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for joint cline, mirror tails
LL.mTails.Joint <- function(pars.JM,x,x.g,n,y)
{p <- Joint.p.mTails(x,x.g,pars.JM)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for spatial cline, left tail
LL.lTails.Joint <- function(pars.JL,x,x.g,n,y)
{p <- Joint.p.lTails(x,x.g,pars.JL)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

# Likelihood for spatial cline, right tail
LL.rTails.Joint <- function(pars.JR,x,x.g,n,y)
{p <- Joint.p.rTails(x,x.g,pars.JR)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}


### Fit models

# Fit Barton spatial cline wiht minmax, no tails, spatial transformed [0,1]
fit.noTails.minmax <- function(x.s,n,y)
{res.noTails <- multistart(parmat=pars.A.minmax,
                           fn=LL.noTails.minmax,
                           x=x.s,n=n,y=y,
                           lower=c(-0.1,0,0,0),
                           upper=c(1.1,2.2,1,1),
                           method="L-BFGS-B",
                           control=list(maximize=TRUE))
est.noTails <- round(res.noTails[which.max(res.noTails$value),],2)
est.noTails}

# Fit tailed joint cline, both tails

fit.bTails.Joint <- function(x.s,x.g,n,y)
{res.bTails <- multistart(par=pars.JB,
                          fn=LL.bTails.Joint,
                          x=x.s,x.g=x.g,n=n,y=y,
                          lower=c(-0.1,0,0,0,0,0,-Inf,0),
                          upper=c(1.1,2.2,2.2,1,2.2,1,Inf,Inf),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.bTails <- round(res.bTails[which.max(res.bTails$value),],2)
est.bTails}

# Fit joint cline, mirrored tails
fit.mTails.Joint <- function(x.s,x.g,n,y)
{res.mTails <- multistart(par=pars.JM,
                          fn=LL.mTails.Joint,
                          x=x.s,x.g=x.g,n=n,y=y,
                          lower=c(-0.1,0,0,0,-Inf,0),
                          upper=c(1.1,2.2,2.2,1,Inf,Inf),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.mTails <- round(res.mTails[which.max(res.mTails$value),],2)
est.mTails}

# Fit joint cline, left tail, spatial transformed [0,1]
fit.lTails.Joint <- function(x.s,x.g,n,y)
{res.lTails <- multistart(par=pars.JL,
                          fn=LL.lTails.Joint,
                          x=x.s,x.g=x.g,n=n,y=y,
                          lower=c(-0.1,0,0,0,-Inf,0),
                          upper=c(1.1,2.2,2.2,1,Inf,Inf),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.lTails <- round(res.lTails[which.max(res.lTails$value),],2)
est.lTails}

# Fit joint cline, right tail, spatial transformed [0,1]
fit.rTails.Joint <- function(x.s,x.g,n,y)
{res.rTails <- multistart(par=pars.JR,
                          fn=LL.rTails.Joint,
                          x=x.s,x.g=x.g,n=n,y=y,
                          lower=c(-0.1,0,0,0,-Inf,0),
                          upper=c(1.1,2.2,2.2,1,Inf,Inf),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
est.rTails <- round(res.rTails[which.max(res.rTails$value),],2)
est.rTails}

# fitting all the models
fitTails.Joint <- function(x,x.g,n,y)
{
  x.s <- trans.down(x)
  
  #Data frame
  tailModels <- data.frame(c=rep(NA,5),w=NA,dL=NA,tL=NA,dR=NA,tR=NA,mu=NA,nu=NA,lnL=NA,k=NA)
  rownames(tailModels) <- c("noTails","rTail","lTail","mTails","bTails")
  
  #No tails
  est.noTails <- fitJoint.noTails(pars.J,x.s,x.g,n,y)
  tailModels[1,] <- c(est.noTails[1:2],2.2,1,2.2,1,est.noTails[3:5],2)
  
  #Right tail
  est.rTails <- fit.rTails.Joint(x.s,x.g,n,y)
  tailModels[2,] <- c(est.rTails[1:2],2.2,1,est.rTails[3:4],est.rTails[5:7],4)
  
  #Left tail
  est.lTails <- fit.lTails.Joint(x.s,x.g,n,y)
  tailModels[3,] <- c(est.lTails[1:4],2.2,1,est.lTails[5:7],4)
  
  #Mirror tails
  est.mTails <- fit.mTails.Joint(x.s,x.g,n,y)
  tailModels[4,] <- c(est.mTails[1:2],est.mTails[3:4],est.mTails[3:4],est.mTails[5:7],4)
  
  #Both tails
  est.bTails <- fit.bTails.Joint(x.s,x.g,n,y)
  tailModels[5,] <- c(est.bTails[c(1:9)],6)
  
  #Re-scale parameters
  tailModels$c <- trans.up(tailModels$c,x)
  tailModels$w <- tailModels$w*diff(range(x))
  tailModels$dL <- tailModels$dL*diff(range(x))
  tailModels$dR <- tailModels$dR*diff(range(x))
  
  #Model fit
  tailModels$AICc <- calc_AICc(tailModels$lnL,tailModels$k,length(x))
  tailModels$dAICc <- round(tailModels$AICc-min(tailModels$AICc),2)
  tailModels$RL <- round(exp(-0.5*(tailModels$dAICc)),2)
  tailModels$AICc_w <- round(tailModels$RL/sum(tailModels$RL),2)
  
  #Model average for dAICc < 2
  tailAvg <- colSums(tailModels[,1:6]*tailModels$AICc_w)
  res <- list(tailAvg,tailModels);names(res) <- c("tailAvg","tailModels")
  return(res)
}
