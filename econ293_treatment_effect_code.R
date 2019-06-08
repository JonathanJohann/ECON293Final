

library(tidyverse)
library(grf)
library(glmnet)
library(ggplot2)
library(amlinear)


n = 200000
set.seed(1234)

X_1 = c(rep(0, 0.3*n), rlnorm(0.7*n, meanlog = 0, sdlog =2))
X_2 = rexp(n = n)
X_3 = rnorm(n)

W =1/(1+exp((-1)*(1 + 0.0000001*X_1 + 3*X_2 +4*X_3))) < 0.2 |1/(1+exp((-1)*(1 + 0.00001*X_1 + 3*X_2 +4*X_3))) > 0.95
Y = X_1 + c(rep(1, 0.3*n), rep(4, 0.7*n))*W + rnorm(n,sd=2)

df <- cbind(X_1,X_2,X_3)

# Scale the covariates as done in the ATE tutorial

df <- df %>% scale()

#
p=3
col <- c("X_1","X_2","X_3")
X_ns = do.call(cbind, lapply(1:p, function(col){matrix(splines::ns(df[,col],df=3), n, 3)}))
dim_ns = dim(X_ns)[2]
X_ns = stats::model.matrix(~.*.-1, data.frame(X_ns)) # pairwise interaction (not including squared term for each column)
X_ns_sq = do.call(cbind, lapply(1:dim_ns, function(col){matrix(X_ns[,col]^2)})) # squared term for each column
X_ns = cbind(X_ns, X_ns_sq)

rm(df)
rm(X_ns_sq)
rm(X_1)
rm(X_2)
rm(X_3)
# 
glmnet.fit.propensity = amlinear:::crossfit.cv.glmnet(X_ns, W,family = "binomial")
theta.hat = amlinear:::crossfit.predict(glmnet.fit.propensity)
p_lasso = 1/(1 + exp(-theta.hat))

indices <- 1:n
lasso_clipped <- indices[(p_lasso>0.05)&(p_lasso<0.95)]
{plot(smooth.spline(p_lasso[lasso_clipped], W[lasso_clipped], df = 4))
  abline(0, 1)
  title("Calibration Plot for Lasso")}


# Fit the propensity model

e.fit <- grf::regression_forest(X=X_ns,Y=W,compute.oob.predictions = TRUE)
e.preds <- predict(e.fit)$predictions
e.fit.fit <- grf::regression_forest(X=matrix(e.preds,ncol=1),Y=W,compute.oob.predictions = TRUE)
e.fit.fit.preds <- predict(e.fit.fit)$predictions
e_hat <- e.fit.fit.preds

# Checking the overlap over buckets
overlap <- function(e,W){
  treated <- e[W==1]
  control <- e[W==0]
  for(i in 1:10){
    ct <- sum((treated>=(i-1)/10) & (treated<(i)/10))
    cc <- sum((control>=(i-1)/10) & (control<(i)/10))
    print(paste((i-1)*10,'% to ',i*10,'% ---- ',ct,' TREATED | ',cc,' CONTROL',sep=""))
  }
}
overlap(e_hat,W)

# Plot the Calibration Plot
rf_clipped <- indices[(e_hat>0.05)&(e_hat<0.95)]
{plot(smooth.spline(e_hat[rf_clipped],W[rf_clipped],df=4))
  abline(0,1)
  title("Calibration Plot for RF")}

# Plot the histogram of fit

prop <- data.frame(W=W[rf_clipped],e_hat=e_hat[rf_clipped])
ggplot(data=prop, aes(x=e_hat, fill=W)) + geom_density(alpha=.3) + ggtitle("Propensity Overlap")

# The larger likelihood, the better
loglik = c(Lasso=mean(W[lasso_clipped] * log(p_lasso[lasso_clipped]) + (1 - W[lasso_clipped]) * log(1 - p_lasso[lasso_clipped])),
           RF=mean(W[rf_clipped] * log(e_hat[rf_clipped]) + (1 - W[rf_clipped]) * log(1 - e_hat[rf_clipped])))

# Ideally, give some sort of reasoning as to why you prefer AIPW. Reasonably, it's because
# of the orthogonal moments argument that allows better performance but make sure you 
# can articulate that well!

# Estimate the ATE using AIPW and Causal Forest
e_hat <- read.csv("prf.csv")
e_hat <- e_hat$p
cf <- causal_forest(X=X_ns,Y=Y,W=W,W.hat=e_hat)

test_calibration(cf)

ate_cf_aipw = average_treatment_effect(cf,subset=rf_clipped,target.sample = "overlap")#e.fit,subset=indices)
tauhat_rf_aipw = c(ATE=ate_cf_aipw["estimate"],
                   lower_ci=ate_cf_aipw["estimate"] - 1.96 * ate_cf_aipw["std.err"],
                   upper_ci=ate_cf_aipw["estimate"] + 1.96 * ate_cf_aipw["std.err"])
tauhat_rf_aipw

# Estimate the ATE using AIPW Lasso


Xmod.int <- X_ns
Xmod.for.lasso = cbind(W, Xmod.int, (2 * W - 1) * Xmod.int)
glmnet.fit.outcome = amlinear:::crossfit.cv.glmnet(Xmod.for.lasso, Y,
                                                   penalty.factor = c(0, rep(1, ncol(Xmod.for.lasso) - 1)))
lasso.yhat.control = amlinear:::crossfit.predict(glmnet.fit.outcome,
                                                 cbind(0, Xmod.int, -Xmod.int))
lasso.yhat.treated = amlinear:::crossfit.predict(glmnet.fit.outcome,
                                                 cbind(1, Xmod.int, Xmod.int))

mean(lasso.yhat.treated - lasso.yhat.control)


indices = rf_clipped
G2 = lasso.yhat.treated[indices] - lasso.yhat.control[indices] +
  W[indices] / e_hat[indices] * (Y[indices] - lasso.yhat.treated[indices]) -
  (1 - W[indices]) / (1 - e_hat[indices]) * (Y[indices] - lasso.yhat.control[indices])#


tau.hat = mean(G2)
se.hat = sqrt(var(G2) / length(G2))
tauhat_lasso_aipw = c(ATE=tau.hat,
                      lower_ci=tau.hat-1.96*se.hat,
                      upper_ci=tau.hat+1.96*se.hat)
tauhat_lasso_aipw

# Estimate the ATE using AIPW Lasso with Approx Residual Debalancing
balancing.weights = amlinear::balance_minimax(Xmod.int[rf_clipped,], W[rf_clipped], zeta = 0.5)
G.balance = lasso.yhat.treated[rf_clipped] - lasso.yhat.control[rf_clipped] +
  balancing.weights * (Y[rf_clipped] - W[rf_clipped] * lasso.yhat.treated[rf_clipped]
                       - (1 - W[rf_clipped]) * lasso.yhat.control[rf_clipped])
tau.hat = mean(G.balance)
se.hat = sqrt(var(G.balance) / length(G.balance))
tauhat_lasso_balance = c(ATE=tau.hat,
                         lower_ci=tau.hat-1.96*se.hat,
                         upper_ci=tau.hat+1.96*se.hat)
tauhat_lasso_balance



# Omnibus test for ATE
mx <- read.csv("yhattest_preds.csv")$x
e_hat <- read.csv("prf.csv")$p
clipped <- which((e_hat<0.95)&(e_hat>0.05))
Y_tilde <- Y- mx
W_tilde <- W-e_hat
Y_tilde <- Y_tilde[clipped]
W_tilde <- W_tilde[clipped]
# 3.189 AIPW-Lasso
# 3.060 Causal Forest
# 3.120 Approx Res
# 3.069 Weighted CF

ate_est <- 3.12
aW <- ate_est * W_tilde
summary(lm(Y_tilde ~ aW -1))

# =====================================================
#
# Heterogeneous treatment effects
#
# =====================================================

# We are going to assume unconfoundedness. Briefly talk about the Instrumental Variable assumption
# and how we cannot entertain this assumption due to lack of domain knowledge. But it's okay
# because we also know the ground truth DGP and so we should be fine

# Ideally, we want to try out T-Learner, S-Learner, X-Learner, and Causal Forest (w/ R Learner)
# Then we want to check for validation of heterogeneous treatment effects followed by
# a comparison using the R-Loss


# Comparing across quartiles
preds <- predict(cf)$predictions
#preds <- preds[which((e_hat<0.95)&(e_hat>0.05))]

summ <- summary(preds[which((e_hat<0.95)&(e_hat>0.05))])
q1_subset <- which(preds<summ[[2]])
q1_subset <- which(q1_subset %in% which((e_hat<0.95)&(e_hat>0.05)))
q2_subset <- which((preds>summ[[2]])&(preds<summ[[3]]))
q2_subset <- which(q2_subset %in% which((e_hat<0.95)&(e_hat>0.05)))
q3_subset <- which((preds>summ[[3]])&(preds<summ[[5]]))
q3_subset <- which(q3_subset %in% which((e_hat<0.95)&(e_hat>0.05)))
q4_subset <- which(preds>summ[[5]])
q4_subset <- which(q4_subset %in% which((e_hat<0.95)&(e_hat>0.05)))


Q1 <- grf::average_treatment_effect(cf,subset=q1_subset,method = "AIPW")
Q2 <- grf::average_treatment_effect(cf,subset=q2_subset,method = "AIPW")
Q3 <- grf::average_treatment_effect(cf,subset=q3_subset,method = "AIPW")
Q4 <- grf::average_treatment_effect(cf,subset=q4_subset,method = "AIPW")



X = X_ns
rm(X_ns)
print("Starting S-Learner...")
fitS <- grf::regression_forest(X=cbind(X,W),Y=Y)
mu1 <- predict(fitS,cbind(X,data.frame(W=rep(1,dim(X)[1]))))
mu0 <- predict(fitS,cbind(X,data.frame(W=rep(0,dim(X)[1]))))
s_preds <- as.matrix(mu1) - as.matrix(mu0)
s_preds <- s_preds[indices]
print("Finished S-Learner...")
print("Starting T-Learner...")
fitT1 <- grf::regression_forest(X=X[W==1,], Y=Y[W==1])
fitT0 <- grf::regression_forest(X=X[W==0,], Y=Y[W==0])
mu1 <- predict(fitT1,X)$predictions 
mu0 <- predict(fitT0,X)$predictions
t_preds <- mu1 - mu0
#t_preds <- t_preds[indices]
write.csv(t_preds,file="Tlearnerpreds.csv")

print("Finished T-Learner...")

print("Starting Causal Forest...")
cf <- grf::causal_forest(X=X, W=W,
                         Y=Y,num.threads=4)
cf.preds <- predict(cf,X)$predictions
cf.preds <- cf.preds[indices]
print("Finished Causal Forest...")

print("Starting X-Learner...")
tf0 = regression_forest(X[W==0,], Y[W==0])
yhat0 = predict(tf0, X[W==1,])$predictions
xf1 = regression_forest(X[W==1,], Y[W==1]-yhat0)
xf.preds.1 = predict(xf1, X)$predictions
# this line ensures we make OOB predictions when appropriate xf.preds.1[W==1] = predict(xf1)$predictions
print("8")
tf1 = regression_forest(X[W==1,], Y[W==1])
yhat1 = predict(tf1, X[W==0,])$predictions
xf0 = regression_forest(X[W==0,], yhat1-Y[W==0])
xf.preds.0 = predict(xf0, X)$predictions
# this line ensures we make OOB predictions when appropriate xf.preds.0[W==0] = predict(xf0)$predictions
print("9")
propf = regression_forest(X, W) 
ehat = predict(propf)$predictions
preds.xf = (1 - ehat) * xf.preds.1 + ehat * xf.preds.0

write.csv(preds.xf,file="XLearnerpreds.csv")

tauhat_s_test <- s_preds 
tauhat_t_test <- t_preds 
tauhat_cf_test <- cf.preds 
tauhat_xl_test <- preds.xf[indices]

#Fitting m_hat and e_hat
Y.forest.test = regression_forest(X = as.matrix(X), Y = Y) 
Y.hat.test = predict(Y.forest.test)$predictions 
W.forest.test = regression_forest(X = as.matrix(X), Y = W) 
W.hat.test = predict(W.forest.test)$predictions

#make sure they're cross fitted

W_tilde <- W- e_hat
Y_tilde <- Y - Y.hat.test
df <- cbind(W_tilde,Y_tilde)



#Randomly shuffle the data

df <- data.frame(df[rf_clipped,])
robinson <- matrix(0,nrow=nrow(df),ncol=1)
yourData<-df[sample(nrow(df)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  trainData <- yourData[-testIndexes, ]
  #Use the test and train data partitions however you desire...
  fit <- lm(Y_tilde ~ W_tilde,data=data.frame(trainData))
  preds <- predict(fit,data.frame(testData))
  robinson[testIndexes,] <- preds
}
robinson <- c(robinson)
write.csv(robinson,file="robinson.csv")
n = length(rf_clipped)
robinson<- fit$coefficients[2][[1]]

tauhat_cf_test <- cf.preds[rf_clipped]
tauhat_s_test <- read.csv("SLearnerpreds.csv")$predictions[rf_clipped]
tauhat_t_test <- read.csv("Tlearnerpreds.csv")$x[rf_clipped]
tauhat_xl_test <- read.csv("XLearnerpreds.csv")$x[rf_clipped]
W.hat.test <- e_hat[rf_clipped]
Y.hat.test2 <- Y.hat.test[rf_clipped]
W2 <- W[rf_clipped]
Y2 <- Y[rf_clipped]
mse_rloss2 <- data.frame(
  S_Learner = (-2/n) * (tauhat_s_test - robinson) * (W2-W.hat.test) * (Y2-Y.hat.test2) + (1/n) * (tauhat_s_test^2 - robinson^2)*(W2-W.hat.test)^2,
  T_Learner = (-2/n) * (tauhat_t_test - robinson) * (W2-W.hat.test) * (Y2-Y.hat.test2) + (1/n) * (tauhat_t_test^2 - robinson^2)*(W2-W.hat.test)^2,
  Causal_Forest = (-2/n) * (tauhat_cf_test - robinson) * (W2-W.hat.test) * (Y2-Y.hat.test2) + (1/n) * (tauhat_cf_test^2 - robinson^2)*(W2-W.hat.test)^2,
  X_Learner = (-2/n) * (tauhat_xl_test - robinson) * (W2-W.hat.test) * (Y2-Y.hat.test2) + (1/n) * (tauhat_xl_test^2 - robinson^2)*(W2-W.hat.test)^2
)

print("Generating mean and sd output...")
mse_rloss_summary2 <- mse_rloss2 %>% summarise_all(.funs=mean) 
means2 <- mse_rloss2 %>% summarise_all(.funs=mean) 
mse_rloss_sd2 <- mse_rloss2 %>% summarise_all(.funs=sd)
sds2 <- mse_rloss_sd2 /sqrt(length(tauhat_s_test))
rloss2 <- cbind(t(means2),t(sds2)) 
colnames(rloss2) <- c("mean","sd") 
print("For R-Loss2")
print(rloss2)
print(rloss2[,1]/rloss2[,2])
saveRDS(rloss2,file=paste("rloss2_",randval,".rds"))








