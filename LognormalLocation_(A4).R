########### Packages to import

require(kernlab)
require(clue)
require(Hotelling)
require(randtoolbox)
require(MASS)
require(LaplacesDemon)
require(GoFKernel)
require(NonNorMvtDist)
require(energy)
require(chi)
require(mvtnorm)
require(truncnorm)
require(transport)
require(ecp)
require(Rcpp)

################### Main function 

mainffinal=function()
{
  altseq=seq(-0.25,-0.01,0.01) # Vector of parameters
  d=2 # controls the dimension
  n=200 # controls the sample size
  alp=0.05 # controls the level of the test
  photel=numeric(length(altseq)) # Usual Hotelling
  penergy=numeric(length(altseq)) # Usual Energy
  pmmd=numeric(length(altseq)) # Usual MMD
  pgrankhotel=numeric(length(altseq)) # Rank Hotelling with Gaussian ERD
  purankhotel=numeric(length(altseq)) # Rank Energy with Gaussian ERD
  pgrankenergy=numeric(length(altseq)) # Rank MMD with Gaussian ERD
  purankenergy=numeric(length(altseq)) # Rank Hotelling with Unif ERD
  pgrankmmd=numeric(length(altseq)) # Rank Energy with Unif ERD
  purankmmd=numeric(length(altseq)) # Rank MMD with Unif ERD
  for(i in 1:length(altseq))
  {
    Sig=diag(d)
    mn=rep(0,d)
    xsam=exp(mvrnorm(n,mu=mn,Sigma=Sig))
    ysam=exp(mvrnorm(n,mu=rep(altseq[i],d),Sigma=Sig))
    photel[i]=as.numeric(hotelling.test(xsam,ysam,perm=TRUE,progBar=FALSE,B=5000)$pval<=alp)
    penergy[i]=as.numeric(eqdist.etest(rbind(xsam,ysam),sizes=c(n,n),R=500)$p.value<=alp)
    pmmd[i]=as.numeric(kmmd(xsam,ysam,asymptotic=TRUE)@AsympH0)
    datam=rbind(xsam,ysam)
    distmat=matrix(0,2*n,2*n)
    hseq=qnorm(halton(2*n,d))
    distmat <- (apply(array(apply(hseq,1,function(x){(x-t(datam))^2}),c(ncol(datam),nrow(datam),nrow(hseq))),2:3,sum))
    assignmentFUN=solve_LSAP(distmat)
    assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
    nch=hseq[assignmentSOL[,2],]
    dm1=nch[(1:n),]
    dm2=nch[-(1:n),]
    pgrankhotel[i]=as.numeric(hotelling.test(dm1,dm2,perm=TRUE,progBar=FALSE,B=5000)$pval<=alp)
    pgrankenergy[i]=as.numeric(eqdist.etest(rbind(dm1,dm2),sizes=c(n,n),R=500)$p.value<=alp)
    pgrankmmd[i]=as.numeric(kmmd(dm1,dm2,asymptotic=TRUE)@AsympH0)
    hsequ=halton(2*n,d)
    distmat <- (apply(array(apply(hsequ,1,function(x){(x-t(datam))^2}),c(ncol(datam),nrow(datam),nrow(hsequ))),2:3,sum))
    assignmentFUN=solve_LSAP(distmat)
    assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
    nch=hsequ[assignmentSOL[,2],]
    dm1=nch[(1:n),]
    dm2=nch[-(1:n),]
    purankhotel[i]=as.numeric(hotelling.test(dm1,dm2,perm=TRUE,progBar=FALSE,B=5000)$pval<=alp)
    purankenergy[i]=as.numeric(eqdist.etest(rbind(dm1,dm2),sizes=c(n,n),R=500)$p.value<=alp)
    purankmmd[i]=as.numeric(kmmd(dm1,dm2,asymptotic=TRUE)@AsympH0)
    print(1000+i)
  }
  oput=rbind(photel,penergy,pmmd,pgrankhotel,pgrankenergy,pgrankmmd,purankhotel,purankenergy,purankmmd,phhg)
  return(oput)
}

############################ Plotting code 

###### Let filename.csv be the data matrix containing the power of the 9 tests above based on the above code
###### Note the above code provides the function for one replication. This was replicated a thousand times parallely.
###### It should have 9 rows and the length of the parameter vector as the number of columns
st=read.csv("filename.csv")[,-1]
xval=seq(0.01,0.2,0.1)
df1=data.frame(x=xval,y=pava(st[1,],dec=T),test = "1")
df2=data.frame(x=xval,y=pava(st[4,],dec=T),test = "2")
df3=data.frame(x=xval,y=pava(st[7,],dec=T),test = "3")
par(bg="grey")
plot(df1$x,df1$y,type="o",pch=11,lwd=2,xlab="Sample size",xaxt="n",ylab="Empirical power",ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,main="d=5",cex.main=2.5)
axis(1,cex.axis=1.5)
lines(df2$x,df2$y,type="o",pch=7,lwd=2,col="red")
lines(df3$x,df3$y,type="o",pch=5,lwd=2,col="blue")
grid(4,4,col="white")
legend("bottomright", bty='n', legend=c("Hotelling","RankGaussHotelling","RankUnifHotelling"),col=c("black","red","blue"), lty=rep(1,6), cex=1.7,lwd=2,pch=c(11,7,5,1),y.intersp=1.5)
