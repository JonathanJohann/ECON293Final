

library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
setwd("/Users/jonathanjohannemann/Documents/ECON293")


rf_clipped <- read.csv("prf.csv")
rf_clipped <- which((rf_clipped$p<0.95)&(rf_clipped$p>0.05))

tauhat_s_test <- read.csv("SLearnerpreds.csv")$predictions[rf_clipped]
tauhat_cf_test <- read.csv("cfpreds.csv")$x[rf_clipped]
tauhat_t_test <- read.csv("Tlearnerpreds.csv")$x[rf_clipped]
tauhat_xl_test <- read.csv("XLearnerpreds.csv")$x[rf_clipped]


s <- ggplot(data=data.frame(CATE=tauhat_s_test),aes(x=CATE,fill=1)) + geom_density(alpha=0.3) + ggtitle("S-Learner HTE")
c<- ggplot(data=data.frame(CATE=tauhat_cf_test),aes(x=CATE,fill=1)) + geom_density(alpha=0.3) + ggtitle("R-Learner (CF) HTE")
t<- ggplot(data=data.frame(CATE=tauhat_t_test),aes(x=CATE,fill=1)) + geom_density(alpha=0.3) + ggtitle("T-Learner HTE")
x<- ggplot(data=data.frame(CATE=tauhat_xl_test),aes(x=CATE,fill=1)) + geom_density(alpha=0.3) + ggtitle("X-Learner HTE")

ggarrange(
  s,c,t,x,ncol=4,
  common.legend = TRUE, legend = "bottom"
)
