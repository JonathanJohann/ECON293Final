

library(tidyverse)
library(grf)
library(glmnet)
library(ggplot2)
library(amlinear)
library(evtree)

n = 200000
set.seed(1234)

X_1 = c(rep(0, 0.3*n), rlnorm(0.7*n, meanlog = 0, sdlog =2))
X_2 = rexp(n = n)
X_3 = rnorm(n)

W =1/(1+exp((-1)*(1 + 0.0000001*X_1 + 3*X_2 +4*X_3))) < 0.2 |1/(1+exp((-1)*(1 + 0.00001*X_1 + 3*X_2 +4*X_3))) > 0.95
Y = X_1 + c(rep(1, 0.3*n), rep(4, 0.7*n))*W + rnorm(n,sd=2)

X <- cbind(X_1,X_2,X_3)

# Scale the covariates as done in the ATE tutorial

X <- X %>% scale()
tau.hat <- read.csv("SLearnerpreds.csv")
prop <- read.csv("prf.csv")
prop.indices <- which((prop$p<0.95)&(prop$p>0.05))
tau.hat <- tau.hat[prop.indices,]


n <- dim(tau.hat)[1]
random_idx <- sample.int(n, size=floor(n/2), replace=F)
cost <- 1.1


Y.hat <- read.csv("yhattest_preds.csv")$x
W.hat <- read.csv("prf.csv")$p
Y.hat <- Y.hat[prop.indices]
W.hat <- W.hat[prop.indices]

Y.hat.train <- Y.hat[random_idx]
W.hat.train <- W.hat[random_idx]
Y.hat.test <- Y.hat[-random_idx]
W.hat.test <- W.hat[-random_idx]
tau.hat.train <- tau.hat[random_idx,"predictions"]
tau.hat.test <- tau.hat[-random_idx,"predictions"]

train.Y <- Y[prop.indices]
train.Y <- train.Y[random_idx]
train.W <- W[prop.indices]
train.W <- train.W[random_idx]

test.Y <- Y[prop.indices]
test.Y <- test.Y[-random_idx]
test.W <- W[prop.indices]
test.W <- test.W[-random_idx]


mu.hat.0 <- Y.hat - W.hat* tau.hat$predictions
mu.hat.1 <- Y.hat + (1 - W.hat) * tau.hat$predictions

# Computing doubly-robust scores
resid <- Y[prop.indices] - W[prop.indices] * mu.hat.1 - (1 - W[prop.indices]) * mu.hat.0
weights <- (W[prop.indices] - W.hat) / (W.hat * (1 - W.hat))

Gamma.hat <- tau.hat$predictions + weights * resid

# We subtract the cost defined in the previous section
Gamma.hat.net <- Gamma.hat - cost


plugin.assignment <- 2*as.numeric(tau.hat > cost) - 1

A.pi <- plugin.assignment*Gamma.hat.net

estimated.benefit <- describe(A.pi)[c("mean", "se")]
estimated.benefit


df_aug <- data.frame(X[prop.indices,])
df_aug <- X[random_idx,]
# Add sign of gamma (denoted Z) and absolute value of gamma (denoted lambda)
df_aug <- cbind(df_aug,data.frame(label=factor(sign(Gamma.hat.train.net))))
df_aug <- cbind(df_aug,data.frame(weights=abs(Gamma.hat.train.net)))

fmla <- as.formula(label ~ X_1 + X_2 + X_3)
opt_policy_tree <- evtree::evtree(formula = fmla, 
                                  data = df_aug,
                                  control=evtree.control(maxdepth=2,
                                                         minbucket=0.025*1000*sum(df_aug$weights),
                                                         minsplit=0.075*1000*sum(df_aug$weights),
                                                         niterations=1000,
                                                         ntrees=100),
                                  weights=round(1000*df_aug$weights))

opt_policy_tree
# Predict optimal assignment

df_train <- data.frame(X[prop.indices,])
df_train <- df_train[random_idx,]
df_test <- data.frame(X[prop.indices,])
df_test <- df_test[-random_idx,]
opt.assignment.train <- as.numeric(as.character(predict(opt_policy_tree, newdata = df_train, type="response")))
opt.assignment.test <- as.numeric(as.character(predict(opt_policy_tree, newdata = df_test, type="response")))

# Calculate value over random policy
opt.A.hat <- data.frame(rbind(
  opt.perf.train = describe(opt.assignment.train*Gamma.hat.train.net)[c("mean", "se")],
  opt.perf.test = describe(opt.assignment.test*Gamma.hat.test.net)[c("mean", "se")]))
opt.A.hat
