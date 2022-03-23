#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("visreg")
#install.packages("emmeans")
#install.packages("sjPlot")
#install.packages("MuMIn")
#install.packages("pbkrtest")
#install.packages("broom.mixed")
#install.packages("gtools")

require(lme4)
require(lmerTest)
require(visreg)
require(emmeans)
require(MuMIn)
require(sjPlot)
require(pbkrtest)
require(broom.mixed)
require(dplyr)

read.csv("/Users/matt/Documents/MATLAB/KRSP/R/RITable.csv")->n
n$season<-as.factor(n$season)
n$mast<-as.factor(n$is_mast)
n$sex<-as.factor(n$sex)
n$preg<-as.factor(n$is_preg)

# remove Winter and add mast interaction
# n_noWinter<-subset(n,n$season!=1)
# see investigate sex correlations in run_this.m
summary(test<-lmer(qb~season*mast+sex*season+age+I(age^2)+(1|squirrel_id),n))
test %>% tidy %>% mutate(signif = stars.pval(p.value))

visreg(test,"sex")
visreg(test,"sex",by="mast")
visreg(test,"sex",by="season")
visreg(test,"age")
visreg(test,"age",trans=exp,partial=TRUE)
visreg(test,"age",by='season',trans=exp,partial=TRUE)
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
emmeans(test, list(pairwise ~ season*sex), adjust = "tukey")
# emmeans(test, list(pairwise ~ mast*sex), adjust = "tukey") # can't do, unless no winter

# see if the results change for QB day or night
summary(test<-lmer(qb_day~season*mast+age+I(age^2)+(1|squirrel_id),n))
test %>% tidy %>% mutate(signif = stars.pval(p.value))
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
visreg(test,"age",by='season')

summary(test<-lmer(qb_night~season*mast+age+I(age^2)+(1|squirrel_id),n))
test %>% tidy %>% mutate(signif = stars.pval(p.value))
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
visreg(test,"age",by='season',trans=exp,partial=TRUE)

n_noMast<-subset(n,n$mast!=1)
n_onlyMast<-subset(n,n$mast!=0)

# this is the pairwise for all seasons irregardless of mast and sex
summary(test<-lmer(qb~season+sex+age+I(age^2)+(1|squirrel_id),n))
emmeans(test, list(pairwise ~ season), adjust = "tukey")
visreg(test,"age",by='season')
visreg(test,'season');

summary(test<-lmer(qb~season+sex+age+I(age^2)+(1|squirrel_id),n_noMast))
emmeans(test, list(pairwise ~ season), adjust = "tukey")
visreg(test,'season');

summary(test<-lmer(qb~season+sex+age+I(age^2)+(1|squirrel_id),n_onlyMast))
emmeans(test, list(pairwise ~ season), adjust = "tukey")
visreg(test,'season');

# does QB depend on trapping incidence?
summary(test<-lmer(qb~season+traps_rec+(1|squirrel_id),n))
visreg(test,"traps_rec",by="season")
summary(test<-lmer(qb_day~traps_rec+(1|squirrel_id),n))
visreg(test,"traps_rec");
summary(test<-lmer(qb_night~traps_rec+(1|squirrel_id),n))
visreg(test,"traps_rec");
# is trapping incidence sex-dependent?
summary(test<-lmer(traps_rec~season*mast+sex+(1|squirrel_id),n))
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")

summary(test<-lmer(qb_nest~season*mast+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"sex",by="season")
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

summary(test<-lmer(trans_per~season*mast+sex+age+I(age^2)+(1|squirrel_id),n))
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
# remove mast
summary(test<-lmer(trans_per~season+sex+age+I(age^2)+(1|squirrel_id),n))
emmeans(test, list(pairwise ~ season), adjust = "tukey")


visreg(test,'season')
visreg(test,"season",by="mast")
visreg(test,"mast",by="season")
visreg(test,"traps_rec",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

n$trans_all<-n$trans_day+n$trans_night
summary(test<-lmer(trans_all~season+(1|squirrel_id),n))
visreg(test,"season")
emmeans(test, list(pairwise ~ season), adjust = "tukey")

summary(test<-lmer(trans_all~season*mast+(1|squirrel_id),n))
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

# removing day/night trans, too much
summary(test<-lmer(trans_day~season*mast+sex+age+I(age^2)++(1|squirrel_id),n))
visreg(test,"trans_day",by="mast")
visreg(test,"trans_day",by="season")
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")
summary(test<-lmer(trans_night~season*mast+sex+age+I(age^2)+(1|squirrel_id),n))
visreg(test,"trans_night",by="mast")
visreg(test,"trans_night",by="season")
visreg(test,"sex",by="season")
visreg(test,"sex",by="mast")

# is RI associated with longevity in the same test?
summary(test<-lmer(RI_odba~season*mast+sex+age+I(age^2)+longevity+(1|squirrel_id),n))
visreg(test,"mast",by="season")
visreg(test,"mast",by="season")
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")

#n_noMast<-subset(n,n$mast!=1)
# !! currently removing seasons 1,2 for midden_cones in MATLAB, no need to do it here
#n_noMast<-subset(n,n$season!=1)
#n_noMast<-subset(n_noMast,n$season!=2)
#n_noMast<-subset(n_noMast,n$season!=3)
summary(test<-lmer(midden_cones~qb+(1|squirrel_id),n)) # age has no effect
visreg(test,"qb")
summary(test<-lmer(midden_cones~RI_odba+(1|squirrel_id),n))
visreg(test,"RI_odba")

summary(test<-lmer(midden_cones_diff~qb+(1|squirrel_id),n))
visreg(test,"qb")
tab_model(test)
summary(test<-lmer(midden_cones_diff~RI_odba+(1|squirrel_id),n))
visreg(test,"RI_odba")
tab_model(test)

summary(test<-lmer(midden_cones_diff~RI_odba+qb+(1|squirrel_id),n))
visreg(test,"RI_odba")
visreg(test,"qb")

## PREGNANT
n_noMast<-subset(n,n$mast!=0)
summary(test<-lmer(qb~season*preg+(1|squirrel_id),n_noMast))
visreg(test,"preg",by='season')
summary(test<-lmer(RI_odba~season*preg+(1|squirrel_id),n))
emmeans(test, list(pairwise ~ preg*season), adjust = "tukey")
visreg(test,"preg",by='season')

n_preg<-subset(n,n$is_preg==1)
n_preg<-subset(n,!is.nan(n$litter_size))
summary(test<-lmer(qb~mast*season+litter_size+(1|squirrel_id),n_preg))
summary(test<-lmer(RI_odba~mast*season+litter_size+(1|squirrel_id),n_preg))

## GROWTH
read.csv("/Users/matt/Documents/MATLAB/KRSP/R/GrowthTable.csv")->m
m$season<-as.factor(m$season)
m$mast<-as.factor(m$is_mast)
summary(test<-lmer(growth~season*mast+RI+qb+(1|squirrel_id),m))
visreg(test,"season")
visreg(test,"RI",'season')
visreg(test,"qb",'season')











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

