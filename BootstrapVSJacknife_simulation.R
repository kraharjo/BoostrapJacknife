#########################################################
#         Stat 454 Final Project Simulation             #
#   							#
#  Topic: Bootstrap and Jackknife variance estimation   #
#         for ratio and regression estimations          #
#							#
#    	       Done by: Kenny Raharjo			#
#########################################################


############################################################
#                     Bootstrap Sampling                   #
############################################################

#---------------------------------------------------------------------------
# We want to use SRSWOR or PPS to select the sample(n) from population(N). -
# (xi,yi), the observed values for unit i, will be selected as a pair.     -
#---------------------------------------------------------------------------

nsim = 1000  # number of  simulations
N = 1000     # population size
n = 100      # sample size

#-------------------------------------------------------------------------------------------
# We want to simulate the population using the model (1 + rexp(N)) and y = b0 + b1*x + ei. -
#-------------------------------------------------------------------------------------------

X = 1 + rexp(N)  # simulated population
e = rnorm(N)     # epsilon, iid random variable with mean 0
b0 = 7           # assuming beta_0 = 7 
b1 = 14          # assuming beta_1 = 14
Y = 7 + 10*X + e # corresponding y-values 



###########################################
# Step (1a) Selecting sample using SRSWOR #
###########################################

SRS = sample(N,n) # Selection using SRSWOR (position)
SRSx = X[SRS]     # x-component of sample
SRSy = Y[SRS]     # y-component of sample

########################################
# Step (1b) Selecting sample using PPS #
########################################

syspps=function(x,n){
N=length(x)
U=sample(N,N)
xx=x[U]
z=rep(0,N)
for(i in 1:N) z[i]=n*sum(xx[1:i])/sum(x)
r=runif(1)
s=numeric()
for(i in 1:N){
if(z[i]>=r){
s=c(s,U[i])
r=r+1
         }
           }
return(s[order(s)])
     }

PPS = syspps(X,n) # Selection using PPS (position)
PPSx = X[PPS]     # x-component
PPSy = Y[PPS]     # y-component



###############################################
# Step (2a) Selecting bootstrap samples (SRS) #
###############################################

bootSRS = sample(n,n,replace=T)         # bootstrap sample
bootSRSx = SRSx[sample(n,n,replace=T)]  # x-component
bootSRSy = SRSy[sample(n,n,replace=T)]  # y-component

###############################################
# Step (2b) Selecting bootstrap samples (PPS) #
###############################################

bootPPS = sample(n,n,replace=T)         # bootstrap sample
bootPPSx = PPSx[sample(n,n,replace=T)]  # x-component
bootPPSy = PPSy[sample(n,n,replace=T)]  # y-component



################################################
# Step (3) Linearization for bootstrap miu_reg #
################################################
#                                              -
#-----------------------------------------------
#        Regression estimator                  -
#-----------------------------------------------
#                                              -
# miu_reg = y_ba + B1 (miu_x - x_ba)           -
# y_ba = mean(y)                               -
# B1 = cov(x,y)/var(x)                         -
# miu_x = population mean = mean(X)            -
# x_ba = sample mean = mean(x)                 -
#                                              -
#-----------------------------------------------
#        Ratio estimator                       -
#-----------------------------------------------
#                                              -
# miu_R = y_ba/x_ba*miu_x                      -
# y_ba = mean(y)                               -
# x_ba = mean(x)                               -
# miu_x = mean(X)                              -
#                                              -
#-----------------------------------------------

miu_reg = function(x,y) mean(y) + cov(x,y)/var(x)*(mean(X)-mean(x))

miu_R = function(x,y) mean(y)/mean(x)*mean(X)



##################################
# Step (4a) Bootstrap loop (SRS) #
##################################

#-----------------------------
# Using regression estimator -
#-----------------------------

SRSreg<-numeric(nsim) # Array to store value of bootstrap

for (i in 1:nsim) {
temp <- sample(n,n,replace=T);
SRSreg[i] = miu_reg (SRSx[temp], SRSy[temp])
}

#------------------------
# Using ratio estimator -
#------------------------

SRSR<-numeric(nsim) # Array to store value of bootstrap

for (i in 1:nsim) {
temp2 <- sample(n,n,replace=T);
SRSR[i] = miu_R (SRSx[temp2], SRSy[temp2])
}

##################################
# Step (4b) Bootstrap loop (PPS) #
##################################

#-----------------------------
# Using regression estimator -
#-----------------------------

PPSreg<- numeric(nsim) # Array to store value of bootstrap

for (i in 1:nsim) {
temp3 <- sample(n,n,replace=T);
PPSreg[i] = miu_reg (PPSx[temp3], PPSy[temp3])
}

#------------------------
# Using ratio estimator -
#------------------------

PPSR<- numeric(nsim) # Array to store value of bootstrap

for (i in 1:nsim) {
temp4 <- sample(n,n,replace=T);
PPSR[i] = miu_R (PPSx[temp4], PPSy[temp4])
}



#############################################################
# Step (5a) Variance estimate and confidence interval (SRS) #
#############################################################

#-----------------------------
# Using regression estimator -
#-----------------------------

SRSregmean = mean (SRSreg)
SRSregvar = var(SRSreg) # variance estimate
SRSregCI = c(SRSregmean-1.96*sqrt(SRSregvar), SRSregmean+1.96*sqrt(SRSregvar)) # 95% confidence interval

#------------------------
# Using ratio estimator -
#------------------------

SRSRmean = mean(SRSR)
SRSRvar = var(SRSR) # variance estimate
SRSRCI = c(SRSRmean-1.96*sqrt(SRSRvar), SRSRmean+1.96*sqrt(SRSRvar)) # 95% confidence interval

#############################################################
# Step (5b) Variance estimate and confidence interval (PPS) #
#############################################################


#-----------------------------
# Using regression estimator -
#-----------------------------

PPSregmean = mean(PPSreg)
PPSregvar = var(PPSreg) # variance estimate
PPSregCI = c(PPSregmean-1.96*sqrt(PPSregvar), PPSregmean+1.96*sqrt(PPSregvar)) # 95% confidence interval

#------------------------
# Using ratio estimator -
#------------------------

PPSRmean = mean(PPSR)
PPSRvar = var(PPSR) # variance estimate
PPSRCI = c(PPSRmean-1.96*sqrt(PPSRvar), PPSRmean+1.96*sqrt(PPSRvar)) # 95% confidence interval









############################################################
#                     Jackknife Sampling                   #
############################################################

#---------------------------------------------------------------------------
# Jackknife sampling algorithm is very similar to bootstrap sampling,	   -
# with the only difference being the sub-samples and the looping algorithm -
#---------------------------------------------------------------------------



##################################
# Step (4a) Jackknife loop (SRS) #
##################################

#------------------------------------------------------------------------
# Algorithm of Jackknife loop:						-
# 									-
# Comparing the i-th Jackknife sample and j-th element of that sample,	-
# if i = j, we remove the element from the Jackknife sample,		-
# if i != j, we keep that element.					-
#------------------------------------------------------------------------

JK_SRS_x <- numeric(n-1) # Array for memory of x-component
JK_SRS_y <- numeric(n-1) # Array for memory of y-component

#-----------------------------
# Using regression estimator -
#-----------------------------

miu_SRS_reg_JK <- numeric(n) # Array for memory of miu_reg

for(i in 1:n){
for(j in 1:n){
	if      (j < i) ((JK_SRS_x[j] <- SRSx[j])     && (JK_SRS_y[j] <- SRSy[j]))
	else if (j > i) ((JK_SRS_x[j-1] <- SRSx[j])   && (JK_SRS_y[j-1] <- SRSy[j])) }
miu_SRS_reg_JK[i] <- mean (JK_SRS_y) + cov(JK_SRS_x,JK_SRS_y)/var(JK_SRS_x)*(mean(X)-mean(JK_SRS_x))}

#------------------------
# Using ratio estimator -
#------------------------

miu_SRS_R_JK <- numeric(n) # Array for memory of miu_R

for(i in 1:n){
for(j in 1:n){
	if      (j < i) ((JK_SRS_x[j] <- SRSx[j])     && (JK_SRS_y[j] <- SRSy[j]))
	else if (j > i) ((JK_SRS_x[j-1] <- SRSx[j])   && (JK_SRS_y[j-1] <- SRSy[j])) }
miu_SRS_R_JK[i] <- mean(JK_SRS_y)/mean(JK_SRS_x)*mean(X)}

##################################
# Step (4B) Jackknife loop (PPS) #
##################################

JK_PPS_x <- numeric(n-1) # Array for memory of x-component
JK_PPS_y <- numeric(n-1) # Array for memory of y-component

#-----------------------------
# Using regression estimator -
#-----------------------------

miu_PPS_reg_JK <- numeric(n) # Array for memory of miu_reg

for(i in 1:n){
for(j in 1:n){
	if      (j < i) ((JK_PPS_x[j] <- PPSx[j])     && (JK_PPS_y[j] <- PPSy[j]))
	else if (j > i) ((JK_PPS_x[j-1] <- PPSx[j])   && (JK_PPS_y[j-1] <- PPSy[j])) }
miu_PPS_reg_JK[i] <- mean (JK_PPS_y) + cov(JK_PPS_x,JK_PPS_y)/var(JK_PPS_x)*(mean(X)-mean(JK_PPS_x))}

#------------------------
# Using ratio estimator -
#------------------------

miu_PPS_R_JK <- numeric(n) # Array for memory of miu_R

for(i in 1:n){
for(j in 1:n){
	if      (j < i) ((JK_PPS_x[j] <- PPSx[j])     && (JK_PPS_y[j] <- PPSy[j]))
	else if (j > i) ((JK_PPS_x[j-1] <- PPSx[j])   && (JK_PPS_y[j-1] <- PPSy[j])) }
miu_PPS_R_JK[i] <- mean(JK_PPS_y)/mean(JK_PPS_x)*mean(X)}



#############################################################
# Step (5a) Variance estimate and confidence interval (SRS) #
#############################################################

#-----------------------------
# Using regression estimator -
#-----------------------------

JK_SRS_regmean = mean(miu_SRS_reg_JK) 
JK_SRS_regvar = (nsim-1)^2/nsim*var(miu_SRS_reg_JK) # variance estimate
JK_SRS_reg_CI = c(JK_SRS_regmean-1.96*sqrt(JK_SRS_regvar), JK_SRS_regmean+1.96*sqrt(JK_SRS_regvar)) # 95% confidence interval

#------------------------
# Using ratio estimator -
#------------------------

JK_SRS_Rmean = mean(miu_SRS_R_JK)
JK_SRS_Rvar = (nsim-1)^2/nsim*var(miu_SRS_R_JK) # variance estimate
JK_SRS_R_CI = c(JK_SRS_Rmean-1.96*sqrt(JK_SRS_Rvar), JK_SRS_Rmean+1.96*sqrt(JK_SRS_Rvar)) # 95% confidence interval

#############################################################
# Step (5b) Variance estimate and confidence interval (PPS) #
#############################################################

#-----------------------------
# Using regression estimator -
#-----------------------------

JK_PPS_regmean = mean(miu_PPS_reg_JK)
JK_PPS_regvar = (nsim-1)^2/nsim*var(miu_PPS_reg_JK) # variance estimate
JK_PPS_reg_CI = c(JK_PPS_regmean-1.96*sqrt(JK_PPS_regvar), JK_PPS_regmean+1.96*sqrt(JK_PPS_regvar))

#------------------------
# Using ratio estimator -
#------------------------

JK_PPS_Rmean = mean(miu_PPS_R_JK)
JK_PPS_Rvar = (nsim-1)^2/nsim*var(miu_PPS_R_JK) # variance estimate
JK_PPS_R_CI = c(JK_PPS_Rmean-1.96*sqrt(JK_PPS_Rvar), JK_SRS_Rmean+1.96*sqrt(JK_PPS_Rvar)) # 95% confidence interval

