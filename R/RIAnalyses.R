#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("visreg")
#install.packages("emmeans")
#install.packages("sjPlot")

require(lme4)
require(lmerTest)
require(visreg)
require(emmeans)

read.csv("/Users/matt/Documents/MATLAB/KRSP/R/RITable.csv")->n

n$season<-as.factor(n$season)
n$mast<-as.factor(n$is_mast)

n<-subset(n,n$season!=1)

summary(test<-lmer(RI~season*mast+(1|squirrel_id),n))
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ mast*season), adjust = "tukey")

summary(test<-lmer(RI~longevity*season*mast+(1|squirrel_id),n))

summary(test<-lmer(RI~age+(1|squirrel_id),n))

summary(test<-lmer(RI~grid_cone_index*season+(1|squirrel_id),n))

summary(test<-lmer(RI~midden_cones*season+(1|squirrel_id),n))
