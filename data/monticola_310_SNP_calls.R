#Code is modified from Pyron et al. 2023; Systematic Biology
#To extract the original ~300 SNPs for cline analysis
#And examine the HZAR tailed functions

#legacy hzar and HIest installation
#install.packages("../Source/hzar_0.2-5.tar.gz", repos = NULL, type = "source")
#install.packages("../Source/HIest_2.0.tar.gz", repos = NULL, type = "source")

##########################
#GBS Analysis - monticola#
##########################
library(adegenet);library(maps);library(viridis)
library(LEA);library(conStruct);library(hzar)
library(HIest);library(readxl);library(scales)
set.seed(1)

###Read in STRUCTURE file from ipyrad
a <- read.structure("./seal_in_c90.str",
                    n.ind = 71,
                    n.loc = 7809,
                    onerowperind = FALSE,
                    col.lab = 1,
                    col.pop = 0,
                    col.others = 0,
                    row.marknames = 0,
                    NA.char = -9)

###Read in data
dat <- read.csv("./gbs_71_localities.csv",header=T,row.names=1)
pops <- factor(dat$Species[match(row.names(a$tab),dat$specimen)])
a$pop <- pops
xy <- dat[,10:9][match(row.names(a$tab),dat$specimen),]
rownames(xy) <- dat$specimen[match(row.names(a$tab),dat$specimen)]
a@other$xy <- xy
a.orig <- a


###Dimensions
#Drop multi-allelic loci & singletons
a <- a[loc=-which(a$loc.n.all>2),drop=TRUE]                                                 #Trim to full biallelic SNPs

#Drop invariant loci
inv <- which(apply(a$tab,2,function(x){length(unique(na.omit(x)))}) %in% c(0,1))            #Trim to taxa, ID invariant loci
qq <- colnames(a$tab[,inv])
vv <- strsplit(qq,"[.]")
zz <- unique(unlist(lapply(vv,'[[',1)))
inv.drop <- match(zz,names(a$all.names))
a <- a[loc=-inv.drop, drop = TRUE]


###Missingness
#SNPs
  x <- a$tab
snp.m <- sort(apply(x,2,function(x){length(which(is.na(x)))})/dim(x)[1])                    #SNP missingness
barplot(snp.m,las=2,ylim=c(0,1),axisnames=FALSE,col="blue")
drop.snp <- names(which(snp.m>0.2)); length(drop.snp)                                       #Plot and count >20% missing
abline(h=0.2,lty=2,col="red")

  #Drop highly missing SNPs
  xx <- strsplit(drop.snp,"[.]")
  yy <- unique(unlist(lapply(xx,'[[',1)))
  snp.drop <- match(yy,names(a$all.names))
  a <- a[loc=-snp.drop, drop = TRUE]
  
  
###Individuals
  x <- a$tab
ind.m <- sort(apply(x,1,function(x){length(which(is.na(x)))})/dim(x)[2])                    #Individual missingness
barplot(ind.m,las=2,ylim=c(0,1),cex.names=0.5,col="blue")
drop.ind <- names(which(ind.m>0.5)); drop.ind                                               #Plot and name >50%
abline(h=0.5,lty=2,col="red")                                                               #Keep all monticola?


###Frequencies
b <- makefreq(a)                                                                            #Recode frequencies
b <- b[,seq(1,dim(b)[2],2)]                                                                 #Trim to major alleles
b <- b[,names(which(apply(b,2,function(x){diff(range(na.omit(x)))})==1))]                   #Trim to complete SNPs
colnames(b) <- unlist(lapply(strsplit(colnames(b),"[.]"),'[[',1))                           #Rename to SNP


###Get cline data
jasper <- c(34.469705609398545, -84.42897508781253)

#cline distances
cline.dists <- vector()
  for (i in 1:dim(xy)[1]) {cline.dists[i] <- hzar.map.greatCircleDistance(xy[i,2], xy[i,1], jasper[1], jasper[2], units = "Km", degrees = TRUE)}
names(cline.dists) <- rownames(xy); cline.dists[which(a$pop == "cheaha")] <- -cline.dists[which(a$pop == "cheaha")]

#rescale SNPs to monticola
cheaha.snps <- which(apply(b,2,FUN=function(x){summary(lm(x~cline.dists))$coef[2,1]})<0)#the cheaha SNPs
c <- b;c[,cheaha.snps] <- 1-c[,cheaha.snps]

#Write SNPs
#write.csv(b,"monticola_310_SNP_freqs.csv")
#write.csv(c,"monticola_310_SNP_freqs_rescaled.csv")
