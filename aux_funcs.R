###########################
#AUXILIARY CLINE FUNCTIONS#
###########################

### Data/parameter transformation functions
trans.down <- function(x,x0=NULL)
    {if(is.null(x0)){x0<-x}      
    (x+abs(min(x0)))/diff(range(x0))}
trans.up <- function(x,x0){(x*diff(range(x0)))-abs(min(x0))}

### AICc calculation
calc_AICc <- function(LL, k, n)
  {AICc <- -2 * LL + 2 * k * (1 + 1 / (n - k - 1))
  return(AICc)}


### Starting parameters

#Barton spatial admixture cline, untransformed
pars.A <- data.frame(c=runif(100,-100,100),
                     w=runif(100,0,100))

#Barton spatial locus cline, no tails, transformed [0,1]
pars.N <- data.frame(c=runif(100),w=runif(100,0,2))

#Barton spatial locus cline, right tail, transformed [0,1]
pars.R <- data.frame(c=runif(100),
                     w=runif(100,0,2),
                     d=runif(100,0,2),
                     t=runif(100))

#Barton spatial locus cline, left tail, transformed [0,1]
pars.L <- data.frame(c=runif(100),
                     w=runif(100,0,2),
                     d=runif(100,0,2),
                     t=runif(100))

#Barton spatial locus cline, mirror tails, transformed [0,1]
pars.M <- data.frame(c=runif(100),
                     w=runif(100,0,2),
                     d=runif(100,0,2),
                     t=runif(100))

#Barton spatial locus cline, both tails, transformed [0,1]
pars.B <- data.frame(c=runif(100),
                     w=runif(100,0,2),
                     dL=runif(100,0,2),
                     tL=runif(100),
                     dR=runif(100,0,2),
                     tR=runif(100))

#Fitzpatrick logit-logistic genomic cline, 
pars.G <- data.frame(mu=runif(100,-100,100),
                     nu=runif(100,0,1000))

#Joint spatial/genomic cline, no tails, spatial transformed [0,1]
pars.J <- data.frame(c=runif(1000),
                     w=runif(1000,0,2),
                     mu=runif(1000,-100,100),
                     nu=runif(1000,0,1000))


### Probability densities
#PDF for the Barton spatial cline, no tails
Barton.p.noTails <- function(x,pars){1/(1 + exp(-((x - pars[1]) * 4/pars[2])))}

#PDF for the Barton spatial cline, right tail
Barton.p.rTails <- function(x,pars.R)
  {ifelse(pars.R[3] < x - pars.R[1], 1 - 1/(1 + 
       exp(4 * pars.R[3]/pars.R[2])) * exp(-(pars.R[4]/(1 + exp(-4 * pars.R[3]/pars.R[2]))) * 
       ((x - pars.R[1]) * 4/pars.R[2] - 4 * pars.R[3]/pars.R[2])), 1/(1 + exp(-((x - 
       pars.R[1]) * 4/pars.R[2]))))}

#PDF for the Barton spatial cline, left tail
Barton.p.lTails <- function(x,pars.L)
{ifelse(pars.L[3] < -(x - pars.L[1]), 1/(1 + 
     exp(4 * pars.L[3]/pars.L[2])) * exp(pars.L[4]/(1 + exp(-4 * pars.L[3]/pars.L[2])) * 
     ((x - pars.L[1]) * 4/pars.L[2] + 4 * pars.L[3]/pars.L[2])), 1/(1 + exp(-((x - 
     pars.L[1]) * 4/pars.L[2]))))}

#PDF for the Barton spatial cline, mirror tails
Barton.p.mTails <- function(x,pars.M)
{ifelse(pars.M[3] < -(x - pars.M[1]), 1/(1 + 
     exp(4 * pars.M[3]/pars.M[2])) * exp(pars.M[4]/(1 + exp(-4 * pars.M[3]/pars.M[2])) * 
     ((x - pars.M[1]) * 4/pars.M[2] + 4 * pars.M[3]/pars.M[2])), ifelse(pars.M[3] < 
     x - pars.M[1], 1 - 1/(1 + exp(4 * pars.M[3]/pars.M[2])) * exp(-(pars.M[4]/(1 + 
     exp(-4 * pars.M[3]/pars.M[2]))) * ((x - pars.M[1]) * 4/pars.M[2] - 4 * 
     pars.M[3]/pars.M[2])), 1/(1 + exp(-((x - pars.M[1]) * 4/pars.M[2])))))}

#PDF for the Barton spatial cline, both tails
Barton.p.bTails <- function(x,pars.B)
{ifelse(pars.B[3] < -(x - pars.B[1]), 1/(1 + 
     exp(4 * pars.B[3]/pars.B[2])) * exp(pars.B[4]/(1 + exp(-4 * pars.B[3]/pars.B[2])) * 
     ((x - pars.B[1]) * 4/pars.B[2] + 4 * pars.B[3]/pars.B[2])), ifelse(pars.B[5] < 
     x - pars.B[1], 1 - 1/(1 + exp(4 * pars.B[5]/pars.B[2])) * exp(-(pars.B[6]/(1 + 
     exp(-4 * pars.B[5]/pars.B[2]))) * ((x - pars.B[1]) * 4/pars.B[2] - 4 * 
     pars.B[5]/pars.B[2])), 1/(1 + exp(-((x - pars.B[1]) * 4/pars.B[2])))))}
  
#PDF for the Fitzpatrick logit-logistic genomic cline
Fitz.p <- function(x,pars.G)
 {mu <- pars.G[1]
  nu <- pars.G[2]
  (x^nu)/((x^nu)+((1-x)^nu)*exp(mu))}

#PDF for the joint Barton/Fitzpatrick spatial/genomic cline, no tails
Joint.p.noTails <- function(x.s,x.g,pars.J){Barton.p.noTails(x.s,pars.J[1:2])*Fitz.p(x.g,pars.J[3:4])}


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

#Likelihood for discrete genotypes
likelihood.D <- function(p,y,n){dbinom(y,n,p,log=TRUE)}

#Likelihood for spatial cline, no tails
LL.noTails <- function(pars.N,x,n,y)
{p <- Barton.p.noTails(x,pars.N)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for spatial cline, right tail
LL.rTails <- function(pars.R,x,n,y)
{p <- Barton.p.rTails(x,pars.R)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for spatial cline, left tail
LL.lTails <- function(pars.L,x,n,y)
{p <- Barton.p.lTails(x,pars.L)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for spatial cline, mirror tails
LL.mTails <- function(pars.M,x,n,y)
{p <- Barton.p.mTails(x,pars.M)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Likelihood for spatial cline, both tails
LL.bTails <- function(pars.B,x,n,y)
{p <- Barton.p.bTails(x,pars.B)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Fitzpatrick logit-logistic likelihood
LL.Fitz <- function(pars.G,x,n,y)
{p <- Fitz.p(x,pars.G)
p <- replace(p, p <= 0, .Machine$double.xmin)
p <- replace(p, p >= 1, 1 - .Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}

#Joint spatial/genomic (Barton/Fitzpatrick)
#Cline likelihood for no tails
LL.Joint.noTails <- function(pars.J,x.s,x.g,n,y)
{p <- Joint.p.noTails(x.s,x.g,pars.J)
p <- replace(p,p<=0,.Machine$double.xmin)
p <- replace(p,p>=1,1-.Machine$double.neg.eps)
LL <- likelihood.C(p,y,n)
return(sum(LL))}


### Fit models

# Fit Barton spatial cline, no tails, spatial transformed [0,1]
fit.noTails <- function(x.s,n,y)
  {res.noTails <- multistart(parmat=pars.N,
                             fn=LL.noTails,
                             x=x.s,n=n,y=y,
                             lower=c(-0.1,0),
                             upper=c(1.1,2.2),
                             method="L-BFGS-B",
                             control=list(maximize=TRUE))
  est.noTails <- round(res.noTails[which.max(res.noTails$value),],2)
  est.noTails}

# Fit Barton spatial cline, right tail, spatial transformed [0,1]
fit.rTails <- function(x.s,n,y)
  {res.rTails <- multistart(par=pars.R,
                           fn=LL.rTails,
                           x=x.s,n=n,y=y,
                           lower=c(-0.1,0,0,0),
                           upper=c(1.1,2.2,2.2,1),
                           method="L-BFGS-B",
                           control=list(maximize=TRUE))
  est.rTails <- round(res.rTails[which.max(res.rTails$value),],2)
  est.rTails}

# Fit Barton spatial cline, left tail, spatial transformed [0,1]
fit.lTails <- function(x.s,n,y)
  {res.lTails <- multistart(par=pars.L,
                         fn=LL.lTails,
                         x=x.s,n=n,y=y,
                         lower=c(-0.1,0,0,0),
                         upper=c(1.1,2.2,2.2,1),
                         method="L-BFGS-B",
                         control=list(maximize=TRUE))
  est.lTails <- round(res.lTails[which.max(res.lTails$value),],2)
  est.lTails}

# Fit Barton spatial cline, mirror tails, spatial transformed [0,1]
fit.mTails <- function(x.s,n,y)
  {res.mTails <- multistart(par=pars.M,
                            fn=LL.mTails,
                            x=x.s,n=n,y=y,
                            lower=c(-0.1,0,0,0),
                            upper=c(1.1,2.2,2.2,1),
                            method="L-BFGS-B",
                            control=list(maximize=TRUE))
    est.mTails <- round(res.mTails[which.max(res.mTails$value),],2)
    est.mTails}

# Fit Barton spatial cline, mirror tails, spatial transformed [0,1]
fit.bTails <- function(x.s,n,y)
  {res.bTails <- multistart(par=pars.B,
                          fn=LL.bTails,
                          x=x.s,n=n,y=y,
                          lower=c(-0.1,0,0,0,0,0),
                          upper=c(1.1,2.2,2.2,1,2.2,1),
                          method="L-BFGS-B",
                          control=list(maximize=TRUE))
  est.bTails <- round(res.bTails[which.max(res.bTails$value),],2)
  est.bTails}

# Fit all Barton spatial cline tail models, in original parameter units
fitTails <- function(x,n,y)
{
  #Spatial transform
  x.s <- trans.down(x)
  
  #Data frame
  tailModels <- data.frame(c=rep(NA,5),w=NA,dL=NA,tL=NA,dR=NA,tR=NA,lnL=NA,k=NA)
  rownames(tailModels) <- c("noTails","rTail","lTail","mTails","bTails")
  
  #No tails
  est.noTails <- fit.noTails(x.s,n,y)
  tailModels[1,] <- c(est.noTails[1:2],2.2,1,2.2,1,est.noTails[3],2)
  
  #Right tail
  est.rTails <- fit.rTails(x.s,n,y)
  tailModels[2,] <- c(est.rTails[1:2],2.2,1,est.rTails[3:4],est.rTails[5],4)
  
  #Left tail
  est.lTails <- fit.lTails(x.s,n,y)
  tailModels[3,] <- c(est.lTails[1:4],2.2,1,est.lTails[5],4)
  
  #Mirror tails
  est.mTails <- fit.mTails(x.s,n,y)
  tailModels[4,] <- c(est.mTails[1:2],est.mTails[3:4],est.mTails[3:4],est.mTails[5],4)
  
  #Both tails
  est.bTails <- fit.bTails(x.s,n,y)
  tailModels[5,] <- c(est.bTails[c(1:7)],6)
  
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

#Fitzpatrick logit-logistic fit
fitFitz <- function(x,n,y)
  {res.Fitz <- multistart(par=pars.G,
                        fn=LL.Fitz,
                        x=x,n=n,y=y,
                        lower=c(-Inf,0.0),
                        upper=c(Inf,Inf),
                        method="L-BFGS-B",
                        control=list(maximize=TRUE))
  est.Fitz <- round(res.Fitz[which.max(res.Fitz$value),],2)
  est.Fitz}

#Barton/Fitzpatrick joint fit
# No tails, spatial/genomic
fitJoint.noTails <- function(pars.J,x.s,x.g,n,y)
  {res.Joint <- multistart(par=pars.J,fn=LL.Joint.noTails,
                            x.s=x.s,x.g=x.g,n=n,y=y,
                            lower=c(-0.1,0,-Inf,0),
                            upper=c(1.1,2.2,Inf,Inf),
                            method="L-BFGS-B",
                            control=list(maximize=TRUE))
  est.J <- round(res.Joint[which.max(res.Joint$value),],2)
  est.J$c <- trans.up(est.J$c,dat$x)
  est.J$w <- est.J$w*diff(range(dat$x))
  est.J}


### Add raster image to plot
addImg <- function
(
    obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
    x = NULL, # mid x coordinate for image
    y = NULL, # mid y coordinate for image
    width = NULL, # width of image (in x coordinate units)
    interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}


#calculate great circle distance using Haversine formula
havDist <- function(lat1, lon1, lat2, lon2) 
{#Convert degrees to radians
  lat1_rad <- lat1 * (pi / 180)
  lon1_rad <- lon1 * (pi / 180)
  lat2_rad <- lat2 * (pi / 180)
  lon2_rad <- lon2 * (pi / 180)
  
  #Differences in latitudes and longitudes
  dlat <- lat2_rad - lat1_rad
  dlon <- lon2_rad - lon1_rad
  
  #Haversine formula
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- 6371.0 * c# Radius of the Earth in kilometers
  
  #Distance in kilometers
  return(distance)}

#distance matrix and distances from farthest point
distances <- function(xy)
  {res <- list()
  distances <- data.frame(matrix(NA,nrow=dim(xy)[1],ncol=dim(xy)[1]),row.names = rownames(xy))
  colnames(distances) <- rownames(distances)
  for(i in 1:dim(distances)[1]){for(j in 1:dim(distances)[1]){distances[i,j] <- havDist(xy[i,1],xy[i,2],xy[j,1],xy[j,2])}}
  res$distances <- distances
  res$dist.orig <- cmdscale(distances,k=1)
  return(res)}

#get hi.mle from SNP matrix
hi.mle <- function(snps)
{(rowSums(snps,na.rm=T)*2)/(dim(snps)[2]*2)}

