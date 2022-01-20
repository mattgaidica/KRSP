#Install packages to do linear mixed effects models. lme4 does the model, lmerTest provides p-values, visreg visualizes the results
# get pakcage info > sessionInfo()
install.packages("lme4")
install.packages("lmerTest")
install.packages("visreg")
install.packages("emmeans")
install.packages("sjPlot")
install.packages('BBmisc')

require(lme4)
require(lmerTest)
require(visreg)
require(emmeans)
require(sjPlot)
require(BBmisc)

#read in data
read.csv("/Users/matt/Documents/MATLAB/KRSP/R/nestAsleepOverlap_v3.csv")->n

#look at data
summary(n)

#make sex, season, mast factors
n$season<-as.factor(n$meanSeason)
n$mast<-as.factor(n$is_mast)
n$sex<-as.factor(n$is_female)


#look at data coverage in mast years vs. non-mast years
n_mast<-subset(n,n$mast=="1")
summary(n_mast$season)

n_mast<-subset(n,n$mast=="0")
summary(n_mast$season)

#no season 1 (winter?) data from mast years so can't compare season1 in mast vs. season1 in non-mast

#################################
# 1. Are there sex differences and seasonal differences in time spent sleeping in nest? 
################################
n$qb<-n$in_asleep+n$out_asleep
summary(test<-lmer(qb~sex+season*mast+(1|squirrelId),n))
visreg(test,"sex")
visreg(test,"season",by="mast")
visreg(test,"season",by="sex")
visreg(test,"sex",by="season")
tab_model(test)
emmeans(test, list(pairwise ~ mast*sex), adjust = "tukey")


#This model does just looks a sex differences and changes across season and includes random effect for squirrel ID
summary(test<-lmer(in_asleep~sex+season*mast+(1|squirrelId),n))
visreg(test,"season",by="mast")
#visualize results
visreg(test)
#season 3 = 4 but lot more time sleeping in nest in winter/spring months
# 

summary(test<-lmer(out_asleep ~sex+season*mast+(1|squirrelId),n))
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
visreg(test)

summary(test<-lmer(in_awake ~sex+season+(1|squirrelId),n))
visreg(test)

summary(test<-lmer(out_awake ~sex+season+(1|squirrelId),n))
################################
# 2. does time spent sleeping in nest vary in mast vs. non-mast years? 
################################
#because there are no season=1 data in mast year, let's exclude them. This model includes an interaction between season 
#and mast to see if there are sig diffs in time spent in nest sleeping in mast vs. non-mast years in seasons 2, 3, and 4

#get rid of season=1
n2<-subset(n,n$season!=1)

summary(test<-lmer(in_asleep ~sex+season*mast+(1|squirrelId),n2))
emmeans(test, list(pairwise ~ season*mast), adjust = "tukey")
#females spend significantly less time in nest sleeping than males
#time spent in nest sleeping is lower in seasons3-4 vs. 2 but not-significantly
#sig interactions between mast and season

#visualize results
visreg(test,"mast",by="season")
#this shows you that time spent in nest sleeping in mast years is slightly lower in seasons 3 and 4 vs. time spent in nest sleeping in 
#seasons 3 and 4 in non-mast years but much of the difference here is much more time spent sleeping in nest in season=2 in mast year than non-mast year

visreg(test,"season",by="mast")
#this highlights how in the one mast year, they are in nest sleeping a lot in season 2 (prob because not breeding) but much less so in season 3 and 4

################################
# 3. does time spent AWAKE in nest vary in mast vs. non-mast years? 
################################
#because there are no season=1 data in mast year, let's exclude them. This model includes an interaction between season and mast to 
#see if there are sig diffs in time spent in nest awake in mast vs. non-mast years in seasons 2, 3, and 4

#get rid of season=1
n2<-subset(n,n$season!=1)

summary(test<-lmer(in_awake ~sex+season*mast+(1|squirrelId),n2))
#no sex diffs
#time spent in nest awake is much lower in season 4 than season 2

#visualize results
visreg(test,"mast",by="season")
#this shows you that time spent in nest awake in mast years is slightly lower in seasons 4 vs. time spent in nest awake in season 4 
#in non-mast years whereas time spent awake in nest in mast years and non-mast years is similar in seasons 2 and 3. This suggests 
#that squirrels in mast years spend less time in the nest overall (sleeping or awake) in autumn than in non-mast years

visreg(test,"season",by="mast")


################################
# 4. making proportion of time in nest spent sleeping
################################
n$total_nest<-n$in_awake+n$in_asleep

#prop nest use spent sleeping
summary(n$in_asleep/n$total_nest)
n$prop.sleep<-(n$in_asleep/n$total_nest)

#seasonal variation in prop time spent sleeping - this is just preliminary, should do glmer model with family=binomial
summary(test<-lmer(prop.sleep~sex+season+(1|squirrelId),n))
#no sex differences
#prop time spent in nest sleeping sig higher as you go from winter to autumn


#does prop time spent in nest sleeping vary according to mast and season?
n2<-subset(n,n$season!=1)
summary(test<-lmer(prop.sleep~sex+season*mast+(1|squirrelId),n2))

visreg(test,"season",by="mast")
visreg(test,"mast",by="season")

emmeans(test, list(pairwise ~ mast*season), adjust = "tukey")

#To do
#ideally the model includes random effect for year of data collection too
#really you want to see how relationship between in_awake vs. in_asleep varies across season and mast yeaars where you include 3-way 
#interaction bretween in_asleep season and mast