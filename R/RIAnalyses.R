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

summary(test<-lmer(RI~season*mast+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"sex",by="season")
emmeans(test, list(pairwise ~ mast*season), adjust = "tukey")

summary(test<-lmer(longevity~RI+(1|byear),n))

#summary(test<-lmer(RI~age+(1|squirrel_id),n))

n2<-subset(n,n$season!=2)
n2<-subset(n,n$season!=3)
summary(test<-lmer(RI~midden_cones*mast+(1|squirrel_id),n2))
summary(test<-lmer(RI~grid_cone_index*mast+(1|squirrel_id),n2))
