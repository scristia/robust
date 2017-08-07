##########################################################################################
####Simulate data
##########################################################################################

##Generate covariate data (complete)
  ###start with 1
m1<-5
s1<-1
nx=100
set.seed(seed = 294)
x1<-rnorm(nx,mean=m1,sd=s1)


##Generate response data (include outliers)
    ###define model
mod_y<-function(x){
  2*x
}

mod_y2<-function(x){
  0.25*pi*sin(pi*x)
}
    ###Noise distribution
      alpha<-0.1
      m0<-0
      sd0<-1
      mout<-0
      sdout<-7


noise_ind<-rbinom(nx,1,alpha)
set.seed(29371)
epsilon<-(1-noise_ind)*rnorm(nx,mean=m0,sd=sd0) + noise_ind*rnorm(nx,mean=mout,sd=sdout)


##Define distribution of R - indicator of missingess
mod_r<-function(x){
  0.3*(sin(5*(x+0.2)) )^2
}
rprob<-mod_r(x1)
set.seed(3587)
r<-rbinom(length(rprob),1,rprob) #20% missing

R<-ifelse(r==1,NA,1)

##Calculate data
y<-R*mod_y(x1) + epsilon
hist(y)


dat1<-data.frame(id=1:100,outlier=noise_ind,noise=epsilon,observed=R,y=y)
filter(dat1,outlier==1)



































