y<-c(5,1,5,14,3,19,1,1,4,22)
t<-c(94,16,63,126,5,31,1,1,2,10)
d<-rbind(y,t)

gibbs <- function(n.sims, beta.start, alpha, gamma, delta,
                  y, t, burnin = 0, thin = 1) {
  beta.draws <- c()
  lambda.draws <- matrix(NA, nrow = n.sims, ncol = length(y))
  beta.cur <- beta.start
  lambda.update <- function(alpha, beta, y, t) {
    rgamma(length(y), y + alpha, t + beta)
  }
  beta.update <- function(alpha, gamma, delta, lambda, y) {
    rgamma(1, length(y) * alpha + gamma, delta + sum(lambda))
  }
  for (i in 1:n.sims) {
    lambda.cur <- lambda.update(alpha = alpha, beta = beta.cur,
                                y = y, t = t)
    beta.cur <- beta.update(alpha = alpha, gamma = gamma,
                            delta = delta, lambda = lambda.cur, y = y)
    if (i > burnin & (i - burnin)%%thin == 0) {
      lambda.draws[(i - burnin)/thin, ] <- lambda.cur
      beta.draws[(i - burnin)/thin] <- beta.cur
    }
  }
  return(list(lambda.draws = lambda.draws, beta.draws = beta.draws))
}

posterior <- gibbs(n.sims = 10000, beta.start = 1, alpha = 1.8,
                   gamma = 0.01, delta = 1, y = y, t = t)
cat("Posterior mean of lambda's :","\n")
colMeans(posterior$lambda.draws)

cat("Posterior mean of beta :","\n")

mean(posterior$beta.draws)

qgamma(c(0.025,0.975),1.8,2.4)

################################################
## Dynamic Pricing with Hierarchical Regression
################################################

library(bayesm)
library(MCMCpack)

data(cheese)
head(cheese)

retailer<-levels(cheese$RETAILER) 
nreg<-length(retailer)

regdata<-NULL 
for (i in 1:nreg) { 
  filter <- cheese$RETAILER==retailer[i] 
  y <- log(cheese$VOLUME[filter]) 
  X <- cbind(1,      # intercept placeholder 
             cheese$DISP[filter], 
             log(cheese$PRICE[filter])) 
  beta.ls[i,] <- lm(y~X-1)$coefficients
  
  regdata[[i]] <- list(y=y, X=X) 
}

Data <- list(regdata=regdata) 
sim.size<-21000
burn.in<-1000
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

beta_draw <- out$betadraw[,1 , burn.in:sim.size]
int <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 ,  burn.in:sim.size]
b_dip <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,3 ,  burn.in:sim.size]
b_price <- round(apply(beta_draw,1,mean),2)

beta_retailer <-cbind.data.frame(retailer,int,b_dip,b_price)
head(beta_retailer)


plot(beta.ls[,2],beta_retailer[,3]
     ,xlab = "leat square estimate of coef for diplay activity"
     ,ylab="posterior mean of coef for diplay activity"
     ,pch=20)


## How much drop in volume is expected if the unit price is increased from 5%? How much will be the effect in revenue?

retailer_n0<-1

retailer_summary <- round(apply(cheese[cheese$RETAILER==retailer[retailer_n0],2:4],2,mean,na.rm=T),3)
retailer_summary

current_price<-2.758
current_offer <- c(1
                   ,retailer_summary["DISP"]
                   ,log(current_price))

alternate_price<-2.62
alternate_offer <- c(1
                     ,retailer_summary["DISP"]
                     ,log(alternate_price))

beta_draw <- out$betadraw[retailer_n0, , 201:2000]
beta_mean <-apply(beta_draw,1,mean)
vol_change<-rep(NA,ncol(beta_draw))
for(i in 1:ncol(beta_draw)){
  volume1 <- exp(beta_draw[,i]%*%current_offer)
  volume2 <- exp(beta_draw[,i]%*%alternate_offer)
  vol_change[i]<-((volume2-volume1)/volume1)*100
}
sum_change <-summary(vol_change)
sum_change
hist(sum_change,col="orange")

expected_revenue = retailer_summary["VOLUME"]*retailer_summary["PRICE"]
names(expected_revenue)<-"current expected revenue"

changed_exp_revenue<-rep(NA,length(vol_change))
for(i in 1:length(vol_change))changed_exp_revenue[i] = retailer_summary["VOLUME"]*(1+vol_change[i]/100)*alternate_price

expected_revenue
summary(changed_exp_revenue)
########################################################
model <- MCMChregress(fixed=log(VOLUME)~DISP+log(PRICE)
                      , random=~DISP+log(PRICE)
                      , group="RETAILER"
                      ,data=cheese
                      , burnin=1000, mcmc=1000, thin=1,verbose=1,seed=123, beta.start=0, sigma2.start=1,Vb.start=1, mubeta=0, Vbeta=1.0E6,r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001)


plot(model$mcmc)

# Summary
summary(model$mcmc)

# Predictive posterior mean for each observation
model$Y.pred

# Predicted-Observed
plot(log(cheese$VOLUME),model$Y.pred)
abline(a=0,b=1,col="red",lwd=2)

plot(cheese$VOLUME,exp(model$Y.pred))
abline(a=0,b=1,col="red",lwd=2)


