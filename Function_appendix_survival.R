### Simulated example data ###

options(scipen=999)
# install.packages("survival")
library(survival)
# install.packages("peperr")
library(peperr)

set.seed(123456)
n <- 1000
m0 <- 0.25; m1 <- -0.5; m2 <- 0.5 #; m3 <- 0; m4 <- 0.2
b0 <- -0.3; b1 <- -0.4; b2 <- 0.5 #; b3 <- 0; b4 <- 0.4
lambdaT <- 0.3 # baseline hazard

a <- rbinom(n,1,0.5)
x <- rbinom(n,1,0.4)
# x1 <- rbinom(n,1,0.4)
# x2 <- rbinom(n,1,0.7)
m <- exp(m0+m1*a+m2*x)/(1+exp(m0+m1*a+m2*x))
m <- rbinom(n,1,m)
# m <- exp(m0+m1*a+m2*x1+m3*x2+m4*a*x2)/(1+exp(m0+m1*a+m2*x1+m3*x2+m4*a*x2))
# m <- rbinom(n,1,m)
T <- rweibull(n,shape=1,scale=1/(lambdaT*exp(b0*a+b1*m+b2*x)))
# T <- rweibull(n,shape=1,scale=1/(lambdaT*exp(b0*a+b1*m+b2*x1+b3*x2+b4*a*x2)))
C <- 6
event <- ifelse(T<C,1,0)
time <- ifelse(event==1,T,C)

data <- data.frame(a,x,m,event,time)
# data <- data.frame(a,x1,x2,m,event,time)

dim(data)
names(data)
head(data)

### Specify outcome model including treatment A (a), mediator M (m), and baseline covariates (x):
# With time = survival or censoring time and event = event indicator

outc.model <- coxph(Surv(time,event) ~ a + m + x,data=data)

summary(outc.model)
plot(survfit(outc.model)$time,survfit(outc.model)$surv,type='l')

### Function ###

eff.mediate.survival <- function(model.y,treat,mediator)
{
data.y <- model.frame(model.y)

data.treat1 <- data.y
data.treat1[,treat]<- 1
data.treat0 <- data.y
data.treat0[,treat]<- 0

## Restricted MLE:
S1 <- survfit(model.y,newdata=data.treat1)
S0 <- survfit(model.y,newdata=data.treat0)
time <- S1$time

ind_0 <- ifelse(data.y[,treat]==0,1,0)
ind_1 <- ifelse(data.y[,treat]==1,1,0)

# Propensity scores as percentages:
ps_0 <- 1-mean(data.y[,treat])
ps_1 <- mean(data.y[,treat])

s1m0.rmle <- apply((t(t(S1$surv)*ind_0/ps_0)),1,mean)
s0m0.rmle <- apply((t(t(S0$surv)*ind_0/ps_0)),1,mean)
s1m1.rmle <- apply((t(t(S1$surv)*ind_1/ps_1)),1,mean)

# Calculate standard errors:

data.treat.base <- data.y

for (i in 2:dim(data.y)[2])
{
	if(is.factor(data.treat.base[,i])) {data.treat.base[,i] <- "0"}
	else {data.treat.base[,i] <- 0}
}
S3 <- survfit(model.y,newdata=data.treat.base)

## SE s1m0.rmle:

# Part one (uncertainty due to propensity score):
U1.S1M0 <- apply(t((t(S1$surv)*ind_0/(ps_0^2))),1,mean)

# Part three (uncertainty due to estimating the baseline hazard):
M <- matrix(nrow=length(data.y[,treat]),ncol=length(time))
ND <- matrix(nrow=length(data.y[,treat]),ncol=length(time))

event <- as.numeric(data.y[,1])[(length(data.y[,treat])+1):(length(data.y[,treat])*2)]

pers_eventtime <- as.numeric(data.y[,1])[1:length(data.y[,treat])]

## Indicator I(T >= u) for every person (rows) at every observed event time (columns):
ind_time <- matrix(nrow=length(data.y[,treat]),ncol=length(time))
for (i in 1:length(data.y[,treat]))
{
  ind_time[i,] <- ifelse((pers_eventtime[i] >= time),1,0)
}

X <- model.matrix(model.y)

Beta <- as.vector(coef(model.y))
b_exp <- as.vector(exp(apply(t(X)*Beta,2,sum))) ## exp(B1E + B2M + B3X)

X1 <- X
X1[,treat]<-1

b_exp_e1 <- as.vector(exp(apply(t(X1)*Beta,2,sum))) ## exp(B1 + B2M + B3X)

a1 <- t(X1)
w <- t(X)

for (i in 1:length(data.y[,treat]))
{
  ## N = 1 if an event is observed, stays 1 thereafter:
  N <- ifelse((pers_eventtime[i] <= time)&(event[i]==1),1,0)
  ND[i,] <- N - c(0,N)[1:length(N)]
  ## Cumulative baseline hazard x exp(B1E + B2M + B3M):
  LE <- -log(survfit(model.y,newdata=data.y[i,])$surv)
  ## Not cummulative for dM(t):
  LE <- LE - c(0,LE)[1:length(LE)]
  M[i,] <- ND[i,] - (ind_time[i,]*LE)
  ## Day after event dM(t) should be zero:
  if (sum(1-N)+1<length(N))
  {
    M[i,(sum(1-N)+2):length(LE)]<-0
  }
}

den <- apply(ind_time*b_exp,2,sum)/length(data.y[,treat])
num <- t(M)
dev <- apply(num/den,2,cumsum)

expect_3_10 <- - ((apply(((t(S1$surv)*ind_0/ps_0) * (b_exp_e1)),2,sum))/length(data.y[,treat]))

U3.S1M0 <- dev*expect_3_10

# Part two (uncertainty due to beta's):

## Score residuals or each individuals contribution to the score vector:
ub <- residuals(model.y,type="score")

## Minus the inverse observed information matrix:
inf <- -solve(cov(ub))

timepart1.10 <- ((t((t(S1$surv)*ind_0/ps_0)*b_exp_e1)*((log(S3$surv))[,1]))%*%t(a1))/length(data.y[,treat])

timepart2.10 <- apply((((t(ind_time*b_exp)%*%t(w))*apply(ND,2,sum))/(apply(ind_time*b_exp,2,sum)^2)),2,cumsum)*
                apply((t(S1$surv)*ind_0/ps_0)*b_exp_e1,2,mean)

U2.S1M0 <- (ub %*% inf) %*% t(timepart1.10 + timepart2.10)

# Part four:
U4.S1M0 <- (t(t(S1$surv)*ind_0/ps_0))

# Variance:
var.s1m0.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s1m0.rmle.1234[i] <- var((U1.S1M0[i]*(data.y[,treat]-ps_1)) + U4.S1M0[i,] +
                                U2.S1M0[,i] + U3.S1M0[i,])
}

SD.logit.s1m0.rmle <- (sqrt(var.s1m0.rmle.1234)/(s1m0.rmle*(1-s1m0.rmle)))

LB.10.rmle <- (s1m0.rmle*exp(-1.96*(SD.logit.s1m0.rmle/sqrt(length(data.y[,treat])))))/
              (1-s1m0.rmle+(s1m0.rmle*exp(-1.96*(SD.logit.s1m0.rmle/sqrt(length(data.y[,treat]))))))
UB.10.rmle <- (s1m0.rmle*exp(+1.96*(SD.logit.s1m0.rmle/sqrt(length(data.y[,treat])))))/
              (1-s1m0.rmle+(s1m0.rmle*exp(+1.96*(SD.logit.s1m0.rmle/sqrt(length(data.y[,treat]))))))

## SE s0m0.rmle:

# Part one (uncertainty due to propensity score):
U1.S0M0 <- apply(t((t(S0$surv)*ind_0/(ps_0^2))),1,mean)

# Part three (uncertainty due to estimating the baseline hazard):
X0 <- X
X0[,treat]<-0

b_exp_e0 <- as.vector(exp(apply(t(X0)*Beta,2,sum))) ## exp(B2M + B3X)

a0 <- t(X0)

expect_3_00 <- - ((apply(((t(S0$surv)*ind_0/ps_0) * (b_exp_e0)),2,sum))/length(data.y[,treat]))

U3.S0M0 <- dev*expect_3_00

# Part two (uncertainty due to beta's):

timepart1.00 <- ((t((t(S0$surv)*ind_0/ps_0)*b_exp_e0)*((log(S3$surv))[,1]))%*%t(a0))/length(data.y[,treat])

timepart2.00 <- apply((((t(ind_time*b_exp)%*%t(w))*apply(ND,2,sum))/(apply(ind_time*b_exp,2,sum)^2)),2,cumsum)*
                apply((t(S0$surv)*ind_0/ps_0)*b_exp_e0,2,mean)

U2.S0M0 <- (ub %*% inf) %*% t(timepart1.00 + timepart2.00)

# Part four:
U4.S0M0 <- (t(t(S0$surv)*ind_0/ps_0))

# Variance:
var.s0m0.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s0m0.rmle.1234[i] <- var((U1.S0M0[i]*(data.y[,treat]-ps_1)) + U4.S0M0[i,] +
                                U2.S0M0[,i] + U3.S0M0[i,])
}
SD.logit.s0m0.rmle <- (sqrt(var.s0m0.rmle.1234)/(s0m0.rmle*(1-s0m0.rmle)))

LB.00.rmle <- (s0m0.rmle*exp(-1.96*(SD.logit.s0m0.rmle/sqrt(length(data.y[,treat])))))/
              (1-s0m0.rmle+(s0m0.rmle*exp(-1.96*(SD.logit.s0m0.rmle/sqrt(length(data.y[,treat]))))))
UB.00.rmle <- (s0m0.rmle*exp(+1.96*(SD.logit.s0m0.rmle/sqrt(length(data.y[,treat])))))/
              (1-s0m0.rmle+(s0m0.rmle*exp(+1.96*(SD.logit.s0m0.rmle/sqrt(length(data.y[,treat]))))))

## SE s1m1.rmle:

# Part one (uncertainty due to propensity score):
U1.S1M1 <- - apply(t((t(S1$surv)*ind_1/(ps_1^2))),1,mean)

# Part three (uncertainty due to estimating the baseline hazard):
expect_3_11 <- - ((apply(((t(S1$surv)*ind_1/ps_1) * (b_exp_e1)),2,sum))/length(data.y[,treat]))

U3.S1M1 <- dev*expect_3_11

# Part two (uncertainty due to beta's):
timepart1.11 <- ((t((t(S1$surv)*ind_1/ps_1)*b_exp_e1)*((log(S3$surv))[,1]))%*%t(a1))/length(data.y[,treat])

timepart2.11 <- apply((((t(ind_time*b_exp)%*%t(w))*apply(ND,2,sum))/(apply(ind_time*b_exp,2,sum)^2)),2,cumsum)*
                apply((t(S1$surv)*ind_1/ps_1)*b_exp_e1,2,mean)

U2.S1M1 <- (ub %*% inf) %*% t(timepart1.11 + timepart2.11)

# Part four:
U4.S1M1 <- (t(t(S1$surv)*ind_1/ps_1))

# Variance:
var.s1m1.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s1m1.rmle.1234[i] <- var((U1.S1M1[i]*(data.y[,treat]-ps_1)) + U4.S1M1[i,] +
                           U2.S1M1[,i] + U3.S1M1[i,])
}
SD.logit.s1m1.rmle <- (sqrt(var.s1m1.rmle.1234)/(s1m1.rmle*(1-s1m1.rmle)))

LB.11.rmle <- (s1m1.rmle*exp(-1.96*(SD.logit.s1m1.rmle/sqrt(length(data.y[,treat])))))/
              (1-s1m1.rmle+(s1m1.rmle*exp(-1.96*(SD.logit.s1m1.rmle/sqrt(length(data.y[,treat]))))))
UB.11.rmle <- (s1m1.rmle*exp(+1.96*(SD.logit.s1m1.rmle/sqrt(length(data.y[,treat])))))/
              (1-s1m1.rmle+(s1m1.rmle*exp(+1.96*(SD.logit.s1m1.rmle/sqrt(length(data.y[,treat]))))))


### Locally efficient estimators

## s1m0.le
U4.ztx <- matrix(nrow=length(time),ncol=length(data.y[,treat]))

second.model <- model.matrix(model.y)[,-grep(treat,colnames(model.matrix(model.y)))]
second.model <- second.model[,-grep(mediator,colnames(second.model))]

data.y$second.predictors <- I(second.model)

for (i in 1:length(time))
{
  U4.ztx[i,] <- predict(glm(S1$surv[i,data.y[,treat]==0] ~ second.predictors,
                data=data.y[data.y[,treat]==0,],
                weights=rep(1/ps_0,times=sum(data.y[,treat]==0)),family=binomial),
                newdata=data.y,type="response")
}

s1m0.le <- apply(U4.ztx,1,sum)/length(data.y[,treat])

U1bis.S1M0 <- U4.ztx/ps_0
U4bis.S1M0 <- U4.ztx/ps_0

var.s1m0.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s1m0.le.1234[i] <- var(((U1.S1M0[i] - U1bis.S1M0[i,])*(data.y[,treat]-ps_1)) +
                               U2.S1M0[,i] + U3.S1M0[i,] + U4.S1M0[i,] +
                              (U4bis.S1M0[i,]*(data.y[,treat]-ps_1)))
}

SD.logit.s1m0.le <- (sqrt(var.s1m0.le.1234)/(s1m0.le*(1-s1m0.le)))

LB.10.le <- (s1m0.le*exp(-1.96*(SD.logit.s1m0.le/sqrt(length(data.y[,treat])))))/
            (1-s1m0.le+(s1m0.le*exp(-1.96*(SD.logit.s1m0.le/sqrt(length(data.y[,treat]))))))
UB.10.le <- (s1m0.le*exp(+1.96*(SD.logit.s1m0.le/sqrt(length(data.y[,treat])))))/
            (1-s1m0.le+(s1m0.le*exp(+1.96*(SD.logit.s1m0.le/sqrt(length(data.y[,treat]))))))


## s0m0.le
U4.ztx0 <- matrix(nrow=length(time),ncol=length(data.y[,treat]))

for (i in 1:length(time))
{
  U4.ztx0[i,] <- predict(glm(S0$surv[i,data.y[,treat]==0] ~ second.predictors,
                 data=data.y[data.y[,treat]==0,],
                 weights=rep(1/ps_0,times=sum(data.y[,treat]==0)),family=binomial),
                 newdata=data.y,type="response")
}

s0m0.le <- apply(U4.ztx0,1,sum)/length(data.y[,treat])

U1bis.S0M0 <- U4.ztx0/ps_0
U4bis.S0M0 <- U4.ztx0/ps_0

var.s0m0.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s0m0.le.1234[i] <- var(((U1.S0M0[i] - U1bis.S0M0[i,])*(data.y[,treat]-ps_1)) +
                               U2.S0M0[,i] + U3.S0M0[i,] + U4.S0M0[i,] +
                              (U4bis.S0M0[i,]*(data.y[,treat]-ps_1)))
}

SD.logit.s0m0.le <- (sqrt(var.s0m0.le.1234)/(s0m0.le*(1-s0m0.le)))

LB.00.le <- (s0m0.le*exp(-1.96*(SD.logit.s0m0.le/sqrt(length(data.y[,treat])))))/
            (1-s0m0.le+(s0m0.le*exp(-1.96*(SD.logit.s0m0.le/sqrt(length(data.y[,treat]))))))
UB.00.le <- (s0m0.le*exp(+1.96*(SD.logit.s0m0.le/sqrt(length(data.y[,treat])))))/
            (1-s0m0.le+(s0m0.le*exp(+1.96*(SD.logit.s0m0.le/sqrt(length(data.y[,treat]))))))


## s1m1.le
U4.ztx1 <- matrix(nrow=length(time),ncol=length(data.y[,treat]))

for (i in 1:length(time))
{
  U4.ztx1[i,] <- predict(glm(S1$surv[i,data.y[,treat]==1] ~ second.predictors,
                 data=data.y[data.y[,treat]==1,],
                 weights=rep(1/ps_1,times=sum(data.y[,treat]==1)),family=binomial),
                 newdata=data.y,type="response")
}

s1m1.le <- apply(U4.ztx1,1,sum)/length(data.y[,treat])

U1bis.S1M1 <- -U4.ztx/ps_1
U4bis.S1M1 <- -U4.ztx/ps_1

var.s1m1.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  var.s1m1.le.1234[i] <- var(((U1.S1M1[i] - U1bis.S1M1[i,])*(data.y[,treat]-ps_1)) +
                               U2.S1M1[,i] + U3.S1M1[i,] + U4.S1M1[i,] +
                               (U4bis.S1M1[i,]*(data.y[,treat]-ps_1)))
}

SD.logit.s1m1.le <- (sqrt(var.s1m1.le.1234)/(s1m1.le*(1-s1m1.le)))

LB.11.le <- (s1m1.le*exp(-1.96*(SD.logit.s1m1.le/sqrt(length(data.y[,treat])))))/
            (1-s1m1.le+(s1m1.le*exp(-1.96*(SD.logit.s1m1.le/sqrt(length(data.y[,treat]))))))
UB.11.le <- (s1m1.le*exp(+1.96*(SD.logit.s1m1.le/sqrt(length(data.y[,treat])))))/
            (1-s1m1.le+(s1m1.le*exp(+1.96*(SD.logit.s1m1.le/sqrt(length(data.y[,treat]))))))

## Direct effect:
# RMLE:

cov.S1M0.S0M0.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  cov.S1M0.S0M0.rmle.1234[i] <- cov(((U1.S1M0[i]*(data.y[,treat]-ps_1)) + U4.S1M0[i,] + U2.S1M0[,i] +
                                      U3.S1M0[i,]),((U1.S0M0[i]*(data.y[,treat]-ps_1)) + U4.S0M0[i,] +
                                      U2.S0M0[,i] + U3.S0M0[i,]))
}

var.direct.ratio.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  B <- matrix(c(var.s1m0.rmle.1234[i],cov.S1M0.S0M0.rmle.1234[i],cov.S1M0.S0M0.rmle.1234[i],var.s0m0.rmle.1234[i]),nrow=2,byrow=F)
  A <- c((1/s0m0.rmle[i]),-(s1m0.rmle[i]/(s0m0.rmle[i])^2))
  var.direct.ratio.rmle.1234[i] <- (A%*%B)%*%t(t(A))
}
var.direct.ratio.rmle.1234.n <- var.direct.ratio.rmle.1234/length(data.y[,treat])

direct.ratio.rmle <- s1m0.rmle / s0m0.rmle
LB.dir.ratio.rmle <- direct.ratio.rmle - (1.96*(sqrt(var.direct.ratio.rmle.1234.n)))
UB.dir.ratio.rmle <- direct.ratio.rmle + (1.96*(sqrt(var.direct.ratio.rmle.1234.n)))

# LE:

cov.S1M0.S0M0.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  cov.S1M0.S0M0.le.1234[i] <- cov(((U1.S1M0[i] - U1bis.S1M0[i,])*(data.y[,treat]-ps_1)) +
                                    U2.S1M0[,i] + U3.S1M0[i,] + U4.S1M0[i,] +
                                   (U4bis.S1M0[i,]*(data.y[,treat]-ps_1)),
                                  ((U1.S0M0[i] - U1bis.S0M0[i])*(data.y[,treat]-ps_1)) +
                                    U2.S0M0[,i] + U3.S0M0[i,] + U4.S0M0[i,] +
                                   (U4bis.S0M0[i,]*(data.y[,treat]-ps_1)))
}

var.direct.ratio.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  B <- matrix(c(var.s1m0.le.1234[i],cov.S1M0.S0M0.le.1234[i],cov.S1M0.S0M0.le.1234[i],var.s0m0.le.1234[i]),nrow=2,byrow=F)
  A <- c((1/s0m0.le[i]),-(s1m0.le[i]/(s0m0.le[i])^2))
  var.direct.ratio.le.1234[i] <- (A%*%B)%*%t(t(A))
}
var.direct.ratio.le.1234.n <- var.direct.ratio.le.1234/length(data.y[,treat])

direct.ratio.le <- s1m0.le/s0m0.le
LB.dir.ratio.le <- direct.ratio.le - (1.96*(sqrt(var.direct.ratio.le.1234.n)))
UB.dir.ratio.le <- direct.ratio.le + (1.96*(sqrt(var.direct.ratio.le.1234.n)))

## Indirect effect:

# RMLE:
cov.S1M0.S1M1.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  cov.S1M0.S1M1.rmle.1234[i] <- cov(((U1.S1M0[i]*(data.y[,treat]-ps_1)) + U4.S1M0[i,] + U2.S1M0[,i] +
                                      U3.S1M0[i,]),((U1.S1M1[i]*(data.y[,treat]-ps_1)) + U4.S1M1[i,] +
                                      U2.S1M1[,i] + U3.S1M1[i,]))
}

var.indirect.ratio.rmle.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  B <- matrix(c(var.s1m1.rmle.1234[i],cov.S1M0.S1M1.rmle.1234[i],cov.S1M0.S1M1.rmle.1234[i],var.s1m0.rmle.1234[i]),nrow=2,byrow=F)
  A <- c((1/s1m1.rmle[i]),-(s1m1.rmle[i]/(s1m0.rmle[i])^2))
  var.indirect.ratio.rmle.1234[i] <- (A%*%B)%*%t(t(A))
}
var.indirect.ratio.rmle.1234.n <- var.indirect.ratio.rmle.1234/length(data.y[,treat])

indirect.ratio.rmle <- s1m1.rmle / s1m0.rmle
LB.indir.ratio.rmle <- indirect.ratio.rmle - (1.96*(sqrt(var.indirect.ratio.rmle.1234.n)))
UB.indir.ratio.rmle <- indirect.ratio.rmle + (1.96*(sqrt(var.indirect.ratio.rmle.1234.n)))

# LE:

cov.S1M0.S1M1.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  cov.S1M0.S1M1.le.1234[i] <- cov(((U1.S1M0[i] - U1bis.S1M0[i,])*(data.y[,treat]-ps_1)) +
                                    U2.S1M0[,i] + U3.S1M0[i,] + U4.S1M0[i,] +
                                   (U4bis.S1M0[i,]*(data.y[,treat]-ps_1)),
                                  ((U1.S1M1[i] - U1bis.S1M1[i])*(data.y[,treat]-ps_1)) +
                                    U2.S1M1[,i] + U3.S1M1[i,] + U4.S1M1[i,] +
                                   (U4bis.S1M1[i,]*(data.y[,treat]-ps_1)))
}

var.indirect.ratio.le.1234 <- rep(NA,times=length(time))
for (i in 1:length(time))
{
  B <- matrix(c(var.s1m1.le.1234[i],cov.S1M0.S1M1.le.1234[i],cov.S1M0.S1M1.le.1234[i],var.s1m0.le.1234[i]),nrow=2,byrow=F)
  A <- c((1/s1m0.le[i]),-(s1m1.le[i]/(s1m0.le[i])^2))
  var.indirect.ratio.le.1234[i] <- (A%*%B)%*%t(t(A))
}
var.indirect.ratio.le.1234.n <- var.indirect.ratio.le.1234/length(data.y[,treat])

indirect.ratio.le <- s1m1.le / s1m0.le
LB.indir.ratio.le <- indirect.ratio.le - (1.96*(sqrt(var.indirect.ratio.le.1234.n)))
UB.indir.ratio.le <- indirect.ratio.le + (1.96*(sqrt(var.indirect.ratio.le.1234.n)))

results <- cbind(s1m0.rmle,LB.10.rmle,UB.10.rmle,s1m0.le,LB.10.le,UB.10.le,s0m0.rmle,LB.00.rmle,UB.00.rmle,s0m0.le,
                 LB.00.le,UB.00.le,s1m1.rmle,LB.11.rmle,UB.11.rmle,s1m1.le,LB.11.le,UB.11.le,direct.ratio.rmle,
                 LB.dir.ratio.rmle,UB.dir.ratio.rmle,direct.ratio.le,LB.dir.ratio.le,UB.dir.ratio.le,indirect.ratio.rmle,
                 LB.indir.ratio.rmle,UB.indir.ratio.rmle,indirect.ratio.le,LB.indir.ratio.le,UB.indir.ratio.le,time)
results

}

### Run eff.mediate.survival function and specify the following:
### 1) outcome model
### 2) treatment name
### 3) mediator name

results <- eff.mediate.survival(model.y=outc.model,treat="a",mediator="m")

### plot example (direct.ratio.rmle vs. direct.ratio.le):
par(mar=c(5.1,4.1,0.5,2.1))
plot(c(results[,"time"],results[,"time"],results[,"time"],results[,"time"],results[,"time"],
       results[,"time"]),c(results[,"direct.ratio.rmle"],results[,"LB.dir.ratio.rmle"],
        results[,"UB.dir.ratio.rmle"],results[,"direct.ratio.le"],
        results[,"LB.dir.ratio.le"],results[,"UB.dir.ratio.le"]), main="", ylab="Risk Ratio",
        type='n', xlab="Time (Years)")
lines(results[,"time"],results[,"direct.ratio.rmle"],lwd=2,lty=1)
lines(results[,"time"],results[,"LB.dir.ratio.rmle"],lwd=1.5,lty=2)
lines(results[,"time"],results[,"UB.dir.ratio.rmle"],lwd=1.5,lty=2)
lines(results[,"time"],results[,"direct.ratio.le"],lwd=2,lty=1,col='blue')
lines(results[,"time"],results[,"LB.dir.ratio.le"],lwd=1.5,lty=2,col='blue')
lines(results[,"time"],results[,"UB.dir.ratio.le"],lwd=1.5,lty=2,col='blue')
legend(0.5,1.6,c("RMLE Direct Effect Ratio","LE Direct Effect Ratio"),
      col=c("black","blue"),cex = 1,lwd=2)

### plot example (indirect.ratio.rmle vs. indirect.ratio.le):
par(mar=c(5.1,4.1,0.5,2.1))
plot(c(results[,"time"],results[,"time"],results[,"time"],results[,"time"],results[,"time"],
       results[,"time"]),c(results[,"indirect.ratio.rmle"],results[,"LB.indir.ratio.rmle"],
        results[,"UB.indir.ratio.rmle"],results[,"indirect.ratio.le"],
        results[,"LB.indir.ratio.le"],results[,"UB.indir.ratio.le"]), main="", ylab="Risk Ratio",
        type='n', xlab="Time (Years)")
lines(results[,"time"],results[,"indirect.ratio.rmle"],lwd=2,lty=1)
lines(results[,"time"],results[,"LB.indir.ratio.rmle"],lwd=1.5,lty=2)
lines(results[,"time"],results[,"UB.indir.ratio.rmle"],lwd=1.5,lty=2)
lines(results[,"time"],results[,"indirect.ratio.le"],lwd=2,lty=1,col='blue')
lines(results[,"time"],results[,"LB.indir.ratio.le"],lwd=1.5,lty=2,col='blue')
lines(results[,"time"],results[,"UB.indir.ratio.le"],lwd=1.5,lty=2,col='blue')
legend(0.5,1.3,c("RMLE Indirect Effect Ratio","LE Indirect Effect Ratio"),
      col=c("black","blue"),cex=1,lwd=2)
