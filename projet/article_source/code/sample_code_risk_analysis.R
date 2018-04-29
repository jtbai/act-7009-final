# author: Antoine J.-P. Tixier
# date: Feb-Mar 2015

# please cite:

# MLA
# Tixier, Antoine J‐P., Matthew R. Hallowell, and Balaji Rajagopalan. "Construction safety risk modeling and simulation." Risk analysis 37.10 (2017): 1917-1935.

# BibTeX
#@article{tixier2017construction,
#  title={Construction safety risk modeling and simulation},
#  author={Tixier, Antoine J-P and Hallowell, Matthew R and Rajagopalan, Balaji},
#  journal={Risk analysis},
#  volume={37},
#  number={10},
#  pages={1917--1935},
#  year={2017},
#  publisher={Wiley Online Library}
#}

# ==============================

library(evmix)

# function courtesy of Prof. Arthur Charpentier
gaussian.kernel.copula.surface=function (u,v,n) {
  s=seq(1/(n+1), length=n, by=1/(n+1))
  mat=matrix(NA,nrow = n, ncol = n)
  sur=kde2d(qnorm(u),qnorm(v),n=1000,lims = c(-4, 4, -4, 4))
  su=sur$z
  for (i in 1:n) {
    for (j in 1:n) {
      Xi=round((qnorm(s[i])+4)*1000/8)+1;
      Yj=round((qnorm(s[j])+4)*1000/8)+1
      mat[i,j]<-su[Xi,Yj]/(dnorm(qnorm(s[i]))*dnorm(qnorm(s[j])))
    }
  }
  return(list(x=s,y=s,z=data.matrix(mat)))
}
library(data.table)
# ========== load data ==========

binary = fread("sample_binary_matrix.csv",header=TRUE)
exp.j = fread("sample_exposure_values.csv",header=TRUE)
j.out = fread("sample_outcomes.csv",header=TRUE,stringsAsFactors=FALSE)

# ========== compute historical risk ==========

# get counts of attributes for each level of real and worst potential outcome
# medical case and lost work time are combined, as usual
counts.level = matrix(nrow=length(colnames(binary)),ncol=10)
colnames(counts.level) = c("real pain","real 1st","real mc.lwt","real pd","real fatality","worst pain","worst 1st","worst mc.lwt","worst pd","worst fatality")
rownames(counts.level)=colnames(binary)
number_of_reports = nrow(binary)

for (i in 1:length(colnames(binary))){
    
    real.occurence = j.out[binary[,i]==1,1]
    worst.occurence = j.out[binary[,i]==1,2]
    
    counts.level[i,] = c(sum(real.occurence=="Pain"),
                       sum(real.occurence=="1st Aid"),
                       sum(real.occurence=="Medical Case")+sum(real.occurence=="Lost Work Time"),
                       sum(real.occurence=="Permanent Disalement"),
                       sum(real.occurence=="Fatality"),
                       sum(worst.occurence=="Pain"),
                       sum(worst.occurence=="1st Aid"),
                       sum(worst.occurence=="Medical Case")+sum(worst.occurence=="Lost Work Time"),
                       sum(worst.occurence=="Permanent Disalement"),
                       sum(worst.occurence=="Fatality"))
}

# sanity check
are_count_equals <- function(actual_by_precursor, expected_by_report){
    all(apply(actual_by_precursor,1,sum)==apply(expected_by_report,2,sum))
}

stopifnot(are_count_equals(counts.level[,1:5], binary))

# severity scores
severity=c(12,48,192,1024,26214)

# risk values
risk = matrix(nrow=length(colnames(binary)),ncol=4)
rownames(risk)=colnames(binary)
colnames(risk)=c("real outcome global risk","worst case scenario global risk","real outcome relative risk","worst case scenario relative risk")

# global risk values
for (i in 1:length(colnames(binary))){
    total.weighted.real.count.by.severity = sum(counts.level[i,1:5]*severity)
    total.weighted.worst.count.by.severity = sum(counts.level[i,6:10]*severity)
    
    risk[i,1] = total.weighted.real.count.by.severity/number_of_reports
    risk[i,2] = total.weighted.worst.count.by.severity/number_of_reports
}

# relative risk values
risk[,3] = risk[,1]/exp.j[,2]
risk[,4] = risk[,2]/exp.j[,2]

# round values up
risk = round(risk,2)

real_not_upscaled_for_frequency = as.matrix(binary)%*%risk[,1]
hist(real_not_upscaled_for_frequency)
# relative risks based on real outcomes
real = as.matrix(binary)%*%risk[,3]
hist(real)

real_risk_severity_per_report = as.matrix(binary)%*%risk[,3]
# relative risks based on worst possible outcomes
worst = as.matrix(binary)%*%risk[,4]
worst_risk_severity_per_report = as.matrix(binary)%*%risk[,4]

# ========== univariate smoothed boostrap with variance correction ==========

X = real_risk_severity_per_report
xeval = seq(min(X), max(X)+sd(X), length=length(X))
neval = length(xeval)
x.bar = mean(X)
var.x = var(X)

# bandwith according to Silverman's rule of thumb
bw = density(X)$bw 

nsim = 1e5
X.sim = as.numeric(vector(length=nsim))
varkern = 1

for (i in 1:nsim){
  X.sim[i] = x.bar+(sample(X,1,replace=TRUE)-x.bar+bw*rnorm(1))/sqrt(1+bw^2*varkern/var.x )
  if(X.sim[i]<0) X.sim[i]=0  # a risk cannot be negative
}

boundary_corrected_domain = seq(from=0,to=max(xeval),length.out=100)
# KDE with boundary correction based on Jones 1993 (local linear fitting at the boundary)
b.c.k = dbckden(boundary_corrected_domain,X.sim,bw=bw,bcmethod="simple")
hist(X,prob=TRUE,main='histogram of historical values with KDE of simulated ones')
lines(x=boundary_corrected_domain,y=b.c.k,lty=1)

quantile(X,0.999)
quantile(X.sim,0.999)

# ========== bivariate analysis ==========

plot(real_risk_severity_per_report,worst_risk_severity_per_report,
     xlab="risk based on real outcomes",
     ylab="risk based on worst possible outcomes",
     main="bivariate construction safety risk",
     cex.axis=0.6,cex.main=0.8,cex.lab=0.7)

grid(lwd=2,col="light grey")

# PSEUDO SPACE

# this transforms the two original RVs to RVs having uniform distributions
ranks_real_risk_severity_per_report = rank(real_risk_severity_per_report)/(length(real_risk_severity_per_report)+1)
ranks_worst_risk_severity_per_report = rank(worst_risk_severity_per_report)/(length(worst_risk_severity_per_report)+1)

U = cbind(ranks_real_risk_severity_per_report, ranks_worst_risk_severity_per_report)

plot(U[,1],U[,2],xlab="pseudo risk based on real outcomes",ylab="pseudo risk based on worst possible outcomes",main="bivariate construction safety risk \n pseudo observations",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)
grid(lwd=2,col="light grey")

# empirical Copula density estimate based on transformed Kernels
output = gaussian.kernel.copula.surface(U[,1],U[,2],n=71)
image(output$x,output$y,output$z,col=gray.colors(50,start=1,end=0),xlab="pseudo risk based on real outcomes",ylab="pseudo risk based on worst potential outcomes",main="nonparametric Copula density estimate",cex.lab=0.5,cex.main=0.8,cex.axis=0.5)
contour(output$x,output$y,output$z,add=TRUE,nlevels=50,drawlabels=FALSE,col="dark grey")
points(U[,1],U[,2])
grid(lwd=2)
box()

# bivariate smoothed bootstrap with variance correction for Copula in original space
real.bar = mean(real)
worst.bar = mean(worst)
var.real = var(real)
var.worst = var(worst)
bw.real = density(real)$bw
bw.worst = density(worst)$bw

varkern = 1
seq = seq(from=1,to=length(real),by=1)

nsim = 1e5
biv.risk.sim = matrix(nrow=nsim,ncol=2)

for (i in 1:nsim){
  
  # randomly select a pair
  j = sample(seq,1)
  x = real[j]
  y = worst[j]
  
  # generate a synthetic pair according to the smoothed boostrap scheme and store it
  x.sim = real.bar+(x-real.bar+bw.real*rnorm(1))/sqrt(1+bw.real^2*varkern/var.real)
  y.sim = worst.bar+(y-worst.bar+bw.worst*rnorm(1))/sqrt(1+bw.worst^2*varkern/var.worst)
  
  biv.risk.sim[i,] = c(x.sim,y.sim)
  
  # a risk cannot be negative, these values will be removed later on
  if(any(biv.risk.sim[i,]<0)) {
    biv.risk.sim[i,]=rep(-999,2)
  }
  
}

# remove -999 values
index.remove = which(apply(biv.risk.sim,1,function(x){any(x==-999)})==TRUE)
biv.risk.sim = biv.risk.sim[-index.remove,]

# compare simulated values to original observations

par(mfrow=c(2,1))
par(mar=rep(4,4))

plot(real,worst,xlab="risk based on real outcomes",ylab="risk based on worst possible outcomes",main="bivariate construction safety risk",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)
grid(lwd=2,col="light grey")

plot(biv.risk.sim[,1],biv.risk.sim[,2],xlab="risk based on real outcomes",ylab="risk based on worst possible outcomes",main="simulated values \n n=10^5",cex.axis=0.6,cex.main=0.8,cex.lab=0.7,pch=".",cex=0.6,xaxt="n")
axis(1,at=seq(0,700,by=100),labels=seq(0,700,by=100),cex.axis=0.6)
abline(v=seq(0,700,by=100),col="lightgray",lty = "dotted",lwd=2)
abline(h=seq(0,10e3,by=2e3),col="lightgray",lty = "dotted",lwd=2)


par(mfrow=c(2,1))
par(mar=rep(4,4))

# bivariate quantile estimation

# conditional quantile
quantile(biv.risk.sim[,1],0.95)
quantile(real,0.95) # just to compare

range(biv.risk.sim[,1])
range(real)

quantile(biv.risk.sim[,2],0.95)
quantile(worst,0.95) # just to compare

range(biv.risk.sim[,2])
range(worst)

# conditional quantile
# say we have evidence that:
real.risk = 200

# select values of worst risk corresponding to real.risk meeting this criteria
conditional = biv.risk.sim[(biv.risk.sim[,1]>=195)&(biv.risk.sim[,1]<=205),2]

quantile(conditional)

# so the nice thing here is that provided evidence, we can provide an estimate of how risky things can get (with a confidence band)

