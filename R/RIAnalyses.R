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
n$sex<-as.factor(n$sex)


# remove Winter and add mast interaction
# n_noWinter<-subset(n,n$season!=1)
n_noMast<-subset(n,n$mast!=1)
summary(test<-lmer(qb~season*mast+sex+age+I(age^2)+longevity+(1|squirrel_id),n))
visreg(test,"sex",by="mast")
visreg(test,"sex",by="season")
visreg(test,"age",'season')
visreg(test,"mast",by="season")
visreg(test,"traps_rec",by="season")
visreg(test,"longevity",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

n_noMast<-subset(n,n$mast!=1)
summary(test<-lmer(qb~season+age+I(age^2)+longevity+(1|squirrel_id),n_noMast))
visreg(test,"age",'season')
visreg(test,"longevity",by="season")

summary(test<-lmer(RI_odba~season+age+I(age^2)+longevity+(1|squirrel_id),n_noMast))
visreg(test,"age",'season')
visreg(test,"longevity",by="season")

# does QB depend on trapping incidence?
n_noMast<-subset(n,n$mast!=0)
summary(test<-lmer(qb~season+traps_rec+(1|squirrel_id),n_noMast))
visreg(test,"traps_rec",by="season")
visreg(test,"traps_rec",by="sex")
# is trapping incidence sex-dependent?
summary(test<-lmer(traps_rec~season*mast+sex+(1|squirrel_id),n))
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")

summary(test<-lmer(qb_nest~season*mast+traps_rec+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"sex",by="season")
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

summary(test<-lmer(trans_per~season*mast+traps_rec+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"season",by="mast")
visreg(test,"mast",by="season")
visreg(test,"traps_rec",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

n$trans_all<-n$trans_day+n$trans_night
summary(test<-lmer(trans_all~season*mast+traps_rec+sex+age+I(age^2)+(1|squirrel_id),n))
# since trapping incidence (traps_rec) predicts QB, is it associated with day/night trans specifically?
summary(test<-lmer(trans_day~season*mast+sex+traps_rec+(1|squirrel_id),n))
visreg(test,"trans_day",by="mast")
visreg(test,"trans_day",by="season")
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")
summary(test<-lmer(trans_night~season*mast+sex+traps_rec+(1|squirrel_id),n))
visreg(test,"trans_night",by="mast")
visreg(test,"trans_night",by="season")
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")


summary(test<-lmer(RI_odba~season*mast+traps_rec+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"mast",by="season")
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

#n_noMast<-subset(n,n$mast!=1)
# !! currently removing seasons 1,2 for midden_cones in MATLAB, no need to do it here
#n_noMast<-subset(n,n$season!=1)
#n_noMast<-subset(n_noMast,n$season!=2)
#n_noMast<-subset(n_noMast,n$season!=3)
summary(test<-lmer(midden_cones~qb+(1|squirrel_id),n))
visreg(test,"qb")
summary(test<-lmer(midden_cones~RI_odba+(1|squirrel_id),n))
visreg(test,"RI_odba")

summary(test<-lmer(midden_cones_diff~qb+(1|squirrel_id),n))
visreg(test,"qb")
summary(test<-lmer(midden_cones_diff~RI_odba+(1|squirrel_id),n))
visreg(test,"RI_odba")

summary(test<-lmer(midden_cones_diff~RI_odba+qb+(1|squirrel_id),n))
visreg(test,"RI_odba")



## NOTES
# this doesn't appear to have enough granularity to trust
summary(test<-lmer(qb~grid_cone_index*mast+(1|squirrel_id),n_noMast))
visreg(test,"grid_cone_index",by="mast")

# this includes mast, but those results are NS
#summary(test<-lmer(RI~age+(1|squirrel_id),n))
n_onlyAutumn<-subset(n,n$season!=1)
n_onlyAutumn<-subset(n,n$season!=2)
n_onlyAutumn<-subset(n,n$season!=3)
summary(test<-lmer(RI_odba~midden_cones*mast+(1|squirrel_id),n_onlyAutumn))
summary(test<-lmer(RI_odba~grid_cone_index*mast+(1|squirrel_id),n_onlyAutumn))
visreg(test,"grid_cone_index",by="mast")

