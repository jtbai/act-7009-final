# author: Antoine J.-P. Tixier
# date: Feb-Mar 2015

# refactorer : Jean-Thomas Baillargeon
# date: April 2018

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

par(mfrow = c(1,1))
# function courtesy of Prof. Arthur Charpentier
gaussian.kernel.copula.surface=function (u,v,n) {
    support = seq(1/(n+1), length=n, by=1/(n+1))
    copula_density = matrix(NA, nrow = n, ncol = n)
    number_of_evaluation = 1000
    min_copula_value = -4
    max_copula_value = 4
    copula_range = max_copula_value - min_copula_value 
    square_copula_limit = c(min_copula_value, max_copula_value, min_copula_value, max_copula_value)
    surface = kde2d(qnorm(u), qnorm(v), n=number_of_evaluation, lims=square_copula_limit)
    surface_density =surface$z
  for (i in 1:n) {
    for (j in 1:n) {
      Xi=round((qnorm(support[i])-min_copula_value)*number_of_evaluation/copula_range)+1;
      Yj=round((qnorm(support[j])-min_copula_value)*number_of_evaluation/copula_range)+1
      copula_density[i,j]<-surface_density[Xi,Yj]/(dnorm(qnorm(support[i]))*dnorm(qnorm(support[j])))
    }
  }
  return(list(x=support,y=support,z=data.matrix(copula_density)))
}
# ========== load data ==========

binary = read.csv("sample_binary_matrix.csv",header=TRUE)
probability_of_precursor_on_site = read.csv("sample_exposure_values.csv",header=TRUE)
j.out = read.csv("sample_outcomes.csv",header=TRUE,stringsAsFactors=FALSE)

# ========== compute historical risk ==========

# get counts of attributes for each level of real and worst potential outcome
# medical case and lost work time are combined, as usual
number_of_precursors = length(colnames(binary))
number_of_reports = nrow(binary)

counts.level = matrix(nrow=number_of_precursors,ncol=10)
colnames(counts.level) = c("real pain","real 1st","real mc.lwt","real pd","real fatality","worst pain","worst 1st","worst mc.lwt","worst pd","worst fatality")
rownames(counts.level) = colnames(binary)

for (precursor_index in 1:number_of_precursors){
    
    real.occurence = j.out[binary[,precursor_index]==1,1]
    worst.occurence = j.out[binary[,precursor_index]==1,2]
    
    counts.level[colnames(binary)[precursor_index],] = c(sum(real.occurence=="Pain"),
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

# risk values
for (precursor_index in 1:number_of_precursors){
    total_weighted_real_count_by_severity = sum(counts.level[precursor_index,1:5]*severity)
    total_weighted_worst_count_by_severity = sum(counts.level[precursor_index,6:10]*severity)
    
    risk[precursor_index,1] = total_weighted_real_count_by_severity/number_of_reports
    risk[precursor_index,2] = total_weighted_worst_count_by_severity/number_of_reports
}

# relative risk values
risk[,3] = risk[,1]/probability_of_precursor_on_site[,2]
risk[,4] = risk[,2]/probability_of_precursor_on_site[,2]

risk = round(risk,2)
real_risk_severity_per_report = as.matrix(binary)%*%risk[,3]
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

hist(X,prob=TRUE,main='histograme des valeurs de risque par rapport')
hist(X,prob=TRUE,main='histogram of historical values with KDE of simulated ones')
lines(x=boundary_corrected_domain,y=b.c.k,lty=1)

#Validation of the methodology
quantile(X,0.999)
quantile(X.sim,0.999)


# ========== bivariate analysis ==========

plot(real_risk_severity_per_report,worst_risk_severity_per_report,
     xlab="risk based on real outcomes",
     ylab="risk based on worst possible outcomes",
     main="bivariate construction safety risk",
     cex.axis=0.6,cex.main=0.8,cex.lab=0.7)


plot(real_risk_severity_per_report,worst_risk_severity_per_report,
     xlab="risque basé sur les résultats réels",
     ylab="risque basé sur les résultats en pire cas",
     main="Risque bivarié",
     cex.axis=0.6,cex.main=0.8,cex.lab=0.7)


grid(lwd=2,col="light grey")

# PSEUDO SPACE

# this transforms the two original RVs to RVs having uniform distributions
ranks_real_risk_severity_per_report = rank(real_risk_severity_per_report)/(length(real_risk_severity_per_report)+1)
ranks_worst_risk_severity_per_report = rank(worst_risk_severity_per_report)/(length(worst_risk_severity_per_report)+1)

U = cbind(ranks_real_risk_severity_per_report, ranks_worst_risk_severity_per_report)


plot(U[,1],U[,2],xlab="pseudo risk based on real outcomes",ylab="pseudo risk based on worst possible outcomes",main="bivariate construction safety risk \n pseudo observations",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)
plot(U[,1],U[,2],xlab="pseudo risque basé sur les résultats réels",ylab="pseudo risque basé sur les résultats en pire cas",main="Risque bivarié \n pseudo observations",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)
grid(lwd=2,col="light grey")

# empirical Copula density estimate based on transformed Kernels
copula_density = gaussian.kernel.copula.surface(U[,1],U[,2],n=71)

#Plotly 3d graph
kd = with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
surface = kde2d(qnorm(U[,1]),qnorm(U[,1]),n=1000,lims = c(-4, 4, -4, 4))

#2dimention KDE visualization
plot_ly(x = surface$x, y = surface$y, z = surface$z) %>% add_surface()

#Copula visualization (3D)
plot_ly(x = copula_density$x, y = copula_density$y, z = copula_density$z) %>% add_surface()

#Copula visualization (2D)
image(copula_density$x,copula_density$y,copula_density$z,col=gray.colors(50,start=1,end=0),xlab="pseudo risk based on real outcomes",ylab="pseudo risk based on worst potential outcomes",main="nonparametric Copula density estimate",cex.lab=0.5,cex.main=0.8,cex.axis=0.5)
image(copula_density$x,copula_density$y,copula_density$z,col=gray.colors(50,start=1,end=0),xlab="pseudo risque basé sur les résultats réels",ylab="pseudo risque basé sur les résultats en pire cas",main="Densité estimé de la copule non paramétrique",cex.lab=0.5,cex.main=0.8,cex.axis=0.5)
contour(copula_density$x,copula_density$y,copula_density$z,add=TRUE,nlevels=50,drawlabels=FALSE,col="dark grey")
points(U[,1],U[,2])
grid(lwd=2)
box()

# bivariate smoothed bootstrap with variance correction for Copula in original space
real.bar = mean(real_risk_severity_per_report)
worst.bar = mean(worst_risk_severity_per_report)
var.real = var(real_risk_severity_per_report)
var.worst = var(worst_risk_severity_per_report)
bw.real = density(real_risk_severity_per_report)$bw
bw.worst = density(worst_risk_severity_per_report)$bw

varkern = 1
support_real_risk = seq(from=1,to=length(real_risk_severity_per_report),by=1)

nsim = 1e5
biv_risk_sim = matrix(nrow=nsim, ncol=2)

for (simulation_index in 1:nsim){
  
  # randomly select a pair
  selected_index = sample(support_real_risk,1)
  x = real_risk_severity_per_report[selected_index]
  y = worst_risk_severity_per_report[selected_index]
  
  # generate a synthetic pair according to the smoothed boostrap scheme and store it
  x_sim = real.bar+(x-real.bar+bw.real*rnorm(1))/sqrt(1+bw.real^2*varkern/var.real)
  y_sim = worst.bar+(y-worst.bar+bw.worst*rnorm(1))/sqrt(1+bw.worst^2*varkern/var.worst)
  
  biv_risk_sim[simulation_index,] = c(x_sim, y_sim)
  
  # a risk cannot be negative, these values will be removed later on
  if(any(biv_risk_sim[simulation_index, ]<0)) {
      biv_risk_sim[simulation_index, ]=rep(-999, 2)
  }
}



# this transforms the two original RVs to RVs having uniform distributions
par(mfrow=c(1,1))
index_to_remove = which(apply(biv_risk_sim, 1, function(x){any(x==-999)})==TRUE)
biv_risk_sim = biv_risk_sim[-index_to_remove, ]

x_simulated =biv_risk_sim[, 1]
y_simulated =biv_risk_sim[, 2]

simulated_ranks_real_risk_severity_per_report = rank(x_simulated)/(length(x_simulated)+1)
simulated_ranks_worst_risk_severity_per_report = rank(y_simulated)/(length(y_simulated)+1)

U_simulated = cbind(simulated_ranks_real_risk_severity_per_report, simulated_ranks_worst_risk_severity_per_report)
# compare simulated values to original observations

par(mfrow=c(2,1))
par(mar=rep(4,4))

#plot(real_risk_severity_per_report,worst_risk_severity_per_report,xlab="risk based on real outcomes",ylab="risk based on worst possible outcomes",main="bivariate construction safety risk",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)
plot(real_risk_severity_per_report,worst_risk_severity_per_report,xlab="risque basé sur les résultats réels",ylab="risque basé sur les résultats en pire cas",main="Risque bivarié ",cex.axis=0.6,cex.main=0.8,cex.lab=0.7)

grid(lwd=2,col="light grey")

#plot(biv_risk_sim[,1],biv_risk_sim[,2],xlab="risk based on real outcomes",ylab="risk based on worst possible outcomes",main="simulated values \n n=10^5",cex.axis=0.6,cex.main=0.8,cex.lab=0.7,pch=".",cex=0.6,xaxt="n")
plot(biv_risk_sim[,1],biv_risk_sim[,2],xlab="risque basé sur les résultats réels",ylab="risque basé sur les résultats en pire cas",main="Valeur simulées \n n=10^5",cex.axis=0.6,cex.main=0.8,cex.lab=0.7,pch=".",cex=0.6,xaxt="n")

axis(1,at=seq(0,700,by=100),labels=seq(0,700,by=100),cex.axis=0.6)
abline(v=seq(0,700,by=100),col="lightgray",lty = "dotted",lwd=2)
abline(h=seq(0,10e3,by=2e3),col="lightgray",lty = "dotted",lwd=2)


par(mfrow=c(2,1))
par(mar=rep(4,4))

# bivariate quantile estimation

# conditional quantile
quantile(biv_risk_sim[,1],0.95)
quantile(real_risk_severity_per_report,0.95) # just to compare

range(biv_risk_sim[,1])
range(real_risk_severity_per_report)

quantile(biv_risk_sim[,2],0.95)
quantile(worst_risk_severity_per_report,0.95) # just to compare

range(biv_risk_sim[,2])
range(worst_risk_severity_per_report)

# conditional quantile

# say we have evidence that:
# The employee was welding overhead and the wind shifted, resulting in discomfort in eye.
real_risk_1 = 10 + 1 + 17
# Worker is unloading a ladder from pickup with bad posture.
real_risk_2 = 15 + 49 + 7 + 3

c = 5  
# select values of worst risk corresponding to real.risk meeting this criteria

# Ajusté pour l'exemple de la présentation
conditional_1 = biv_risk_sim[(biv_risk_sim[,1]>=(real_risk_1-c))&(biv_risk_sim[,1]<=(real_risk_1+c)),2]
quantile(conditional_1, 0.8)

conditional_2 = biv_risk_sim[(biv_risk_sim[,1]>=(real_risk_2-c))&(biv_risk_sim[,1]<=(real_risk_2+c)),2]
quantile(conditional_2, 0.8)

# so the nice thing here is that provided evidence, we can provide an estimate of how risky things can get (with a confidence band)

