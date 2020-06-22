rm(list = ls())
library(bayesm)
library(MCMCpack)
data(cheese)
head(cheese)

retailer<-levels(cheese$RETAILER) 
head(retailer)
#cheese_sub <-cheese[cheese$RETAILER==retailer[i],]

cheese$RETAILER<-as.character(cheese$RETAILER)
cheese$city <- NA
cheese$company <- NA

for(i in 1:nrow(cheese)){
  city_comp <- strsplit(cheese$RETAILER[i], split = "-")
  cheese$city[i]<-city_comp[[1]][1]
  cheese$company[i]<-city_comp[[1]][2]
}

city <- unique(cheese$city)
company <- unique(cheese$company)

north_east <-c("PITTSBURGH ","BOSTON ","HARTFORD ")
mid_west <-c("CHICAGO ","CLEVELAND ","ST. LOUIS ","DENVER ")
south <- c("HOUSTON ","DETROIT ","DALLAS/FT. WORTH ")
west <-c("LOS ANGELES ","SAN FRANCISCO ","TAMPA/ST. PETE ")


cheese_sub <-cheese[cheese$city=="PITTSBURGH ",]
unique(cheese_sub$company)
dim(cheese_sub)


cheese_sub <- rbind(cheese_sub
                    ,cheese[cheese$city=="BOSTON ",])
cheese_sub <- rbind(cheese_sub
                    ,cheese[cheese$city=="NEW YORK (NEW) ",])
cheese_sub <- rbind(cheese_sub
                    ,cheese[cheese$city=="ALBANY,NY ",])

cheese_sub <- rbind(cheese_sub,cheese[cheese$city=="CHICAGO ",])
cheese_sub <- rbind(cheese_sub,cheese[cheese$city=="CLEVELAND ",])
cheese_sub <- data.frame(rbind(cheese_sub,cheese[cheese$city=="ST. LOUIS ",]))

retailer<-unique(cheese_sub$RETAILER) 
nreg<-length(retailer)


regdata<-NULL 
for (i in 1:nreg) { 
  filter <- cheese_sub$RETAILER==retailer[i] 
  y <- log(cheese_sub$VOLUME[filter]) 
  X <- cbind(1,      # intercept placeholder 
             cheese_sub$DISP[filter], 
             log(cheese_sub$PRICE[filter])) 
  regdata[[i]] <- list(y=y, X=X) 
}

Z<-matrix(c(rep(1,8),rep(0,6),rep(0,8),rep(1,6)),ncol=2)
colnames(Z)<-c("New England","Mid-West")
nz = ncol(Z)

Z_ext <- cbind.data.frame(retailer,Z)

sim.size<-2000
burn.in<-floor(sim.size*0.1)
Data <- list(regdata=regdata, Z=Z) 
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

### Retailer level
beta_draw <- out$betadraw[,1 , 201:sim.size]
int <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 , 201:sim.size]
b_dip <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,3 , 201:sim.size]
b_price <- round(apply(beta_draw,1,mean),2)

beta_retailer <-cbind.data.frame(retailer,int,b_dip,b_price)
colnames(beta_retailer)<-c("Retailter","Int","display","Price")
head(beta_retailer)

### Geographical level analysis

Deltadraw <- out$Deltadraw
Deltadraw_star <- matrix(Deltadraw[200,],nrow = nz,byrow = F)
Delta_hat<-apply(Deltadraw[201:sim.size,],2,mean)
Delta_hat<-matrix(Delta_hat,nrow = nz,byrow = F)
colnames(Delta_hat)<-c("Int","display","Price")
rownames(Delta_hat)<-colnames(Z)
Delta_hat
#Z%*%Delta_hat


###########
### a) How much surge in volume is expected if 5% discount on the current unit price is offered? 
### b) How much will be the effect in revenue?
### c) How much will be effect on the profit if cost per unit is $2.00?

unit_cost<-2
retailer_n0<-1

retailer_summary <- round(apply(cheese_sub[cheese_sub$RETAILER==retailer[retailer_n0],2:4],2,mean,na.rm=T),3)
retailer_summary

round(apply(cheese_sub[cheese_sub$RETAILER==retailer[retailer_n0],2:4],2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T),3)

current_price<-2.88
current_offer <- c(1
                   ,retailer_summary["DISP"]
                   ,log(current_price))

alternate_price<-2.74
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
hist(vol_change,col="orange",probability = T)
lines(density(vol_change),col="red",lwd=2)

## effect on Revenue and profit
expected_revenue = retailer_summary["VOLUME"]*retailer_summary["PRICE"]
names(expected_revenue)<-"current expected revenue"

expected_profit <-expected_revenue-(retailer_summary["VOLUME"]*unit_cost)

names(expected_profit)<-"current expected profit"


changed_exp_revenue<-profit<-rep(NA,length(vol_change))
for(i in 1:length(vol_change)){
  changed_exp_revenue[i] = retailer_summary["VOLUME"]*(1+vol_change[i]/100)*alternate_price
  cost <-retailer_summary["VOLUME"]*(1+vol_change[i]/100)*unit_cost
  profit[i] <- changed_exp_revenue[i] - cost
}

expected_revenue
summary(changed_exp_revenue)
quantile(changed_exp_revenue,probs = c(0.025,0.975))
hist(changed_exp_revenue, main="",col="orange")
abline(v=expected_revenue,lwd=2,col="red")

## Change in profit

expected_profit
summary(profit)
quantile(profit,probs = c(0.025,0.975))
hist(profit, main="",xlim = c(2700,3150),col="orange")
abline(v=expected_profit,lwd=2,col="red")

### Analysis for `PITTSBURGH - GIANT EAGLE`

unit_cost<-2
retailer_n0<-8

retailer_summary <- round(apply(cheese_sub[cheese_sub$RETAILER==retailer[retailer_n0],2:4],2,mean,na.rm=T),3)
retailer_summary

round(apply(cheese_sub[cheese_sub$RETAILER==retailer[retailer_n0],2:4],2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T),3)

current_price<-2.88
current_offer <- c(1
                   ,retailer_summary["DISP"]
                   ,log(current_price))

alternate_price<-2.74
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
hist(vol_change,col="orange",probability = T)
lines(density(vol_change),col="red",lwd=2)

## effect on Revenue and profit
expected_revenue = retailer_summary["VOLUME"]*retailer_summary["PRICE"]
names(expected_revenue)<-"current expected revenue"

expected_profit <-expected_revenue-(retailer_summary["VOLUME"]*unit_cost)

names(expected_profit)<-"current expected profit"


changed_exp_revenue<-profit<-rep(NA,length(vol_change))
for(i in 1:length(vol_change)){
  changed_exp_revenue[i] = retailer_summary["VOLUME"]*(1+vol_change[i]/100)*alternate_price
  cost <-retailer_summary["VOLUME"]*(1+vol_change[i]/100)*unit_cost
  profit[i] <- changed_exp_revenue[i] - cost
}

expected_revenue
summary(changed_exp_revenue)
quantile(changed_exp_revenue,probs = c(0.025,0.975))
hist(changed_exp_revenue, main="",col="orange",xlim=c(3750,5000))
abline(v=expected_revenue,lwd=2,col="red")

## Change in profit

expected_profit
summary(profit)
quantile(profit,probs = c(0.025,0.975))
hist(profit, main="",xlim = c(1010,1350),col="orange")
abline(v=expected_profit,lwd=2,col="red")

### How the things are going to affect at Geographical Levels

### Analysis for New England

unit_cost<-2
retailer_n0<-1:8

regional_summary <- round(apply(cheese_sub[cheese_sub$RETAILER==retailer[retailer_n0],2:4],2,mean,na.rm=T),3)
regional_summary

current_price<-2.88
current_offer <- c(1
                   ,regional_summary["DISP"]
                   ,log(current_price))

alternate_price<-2.74
alternate_offer <- c(1
                     ,regional_summary["DISP"]
                     ,log(alternate_price))


Deltadraw_NewEng <- out$Deltadraw[(burn.in+1):sim.size,c(1,3,5)]

Deltadraw_NewEng_mean <-apply(Deltadraw_NewEng,2,mean)
vol_change<-rep(NA,nrow(Deltadraw_NewEng))
for(i in 1:nrow(Deltadraw_NewEng)){
  volume1 <- exp(Deltadraw_NewEng[i,]%*%current_offer)
  volume2 <- exp(Deltadraw_NewEng[i,]%*%alternate_offer)
  vol_change[i]<-((volume2-volume1)/volume1)*100
}
sum_change <-summary(vol_change)
sum_change
hist(vol_change,col="orange"
     ,probability = T
     ,xlab="% Change in Volatility"
     ,sub="New England"
     ,main="")
lines(density(vol_change),col="red",lwd=2)

## effect on Revenue and profit
expected_revenue = regional_summary["VOLUME"]*regional_summary["PRICE"]
names(expected_revenue)<-"current expected revenue"

expected_profit <-expected_revenue-(regional_summary["VOLUME"]*unit_cost)
names(expected_profit)<-"current expected profit"


changed_exp_revenue<-profit<-rep(NA,length(vol_change))
for(i in 1:length(vol_change)){
  changed_exp_revenue[i] = regional_summary["VOLUME"]*(1+vol_change[i]/100)*alternate_price
  cost <-regional_summary["VOLUME"]*(1+vol_change[i]/100)*unit_cost
  profit[i] <- changed_exp_revenue[i] - cost
}

expected_revenue
summary(changed_exp_revenue)
quantile(changed_exp_revenue,probs = c(0.025,0.975))
hist(changed_exp_revenue, main="",col="orange",xlim=c(11500,13500))
abline(v=expected_revenue,lwd=2,col="red")

## Change in profit

expected_profit
summary(profit)
quantile(profit,probs = c(0.025,0.975))
hist(profit, main="",xlim = c(3100,4500),col="orange")
abline(v=expected_profit,lwd=2,col="red")


#### Three Layer Hierarchical Bayes Model

Z_ext
BOSTON<-rep(0,nreg)
BOSTON[2:4]<-1
NY<-rep(0,nreg)
NY[5:7]<-1
CHICAGO<-rep(0,nreg)
CHICAGO[9:11]<-1


Z1<-cbind.data.frame(Z,BOSTON,NY,CHICAGO)
Z1<-as.matrix(Z1)
Z1.ext<-cbind.data.frame(retailer,Z1)

New_England <- apply(out$Deltadraw[,seq(1,9,5)],2,mean)
Mid_West <- apply(out$Deltadraw[,seq(2,15,5)],2,mean)
BOSTON<-apply(out$Deltadraw[,seq(3,15,5)],2,mean)
NY<-apply(out$Deltadraw[,seq(4,15,5)],2,mean)
CHICAGO<-apply(out$Deltadraw[,seq(5,15,5)],2,mean)
post_est_coef <- cbind.data.frame(New_England,Mid_West,BOSTON,NY,CHICAGO)
roenames(post_est_coef) <-c("Int","Disp","Price")

###### with city only model
sim.size<-20000
burn.in<-floor(sim.size*0.1)
Data <- list(regdata=regdata, Z=Z1[,3:5]) 
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

str(out)


BOSTON<-apply(out$Deltadraw[,seq(1,9,3)],2,mean)
NY<-apply(out$Deltadraw[,seq(2,9,3)],2,mean)
CHICAGO<-apply(out$Deltadraw[,seq(3,9,3)],2,mean)
post_est_coef <- cbind.data.frame(BOSTON,NY,CHICAGO)
rownames(post_est_coef) <-c("Int","Disp","Price")

?rank

length(unique(rank(Z)))

length(unique(rank(Z1)))