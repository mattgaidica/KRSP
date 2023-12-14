# Squirrel Axy to Behaviour Script

library(data.table)
library(tidyverse)
library(lubridate)
library(zoo)

#Functions that are needed
tempFilterFall= function(axy){
  ## This function filters the temperature recordings for the spurlious events where a faulty temperature is recorded. It then replaces that temperature with the median temperature between the recording preceeding and following it.  I double checked that these faulty temperature never seem to occur in consequative recordings. Temperature sensor only seems to register a temperature every 10 seconds (probably power saving method), so every 10 recordings have the same temperature.  Once this is used, tempC variable is teh temperature that you will want to use for all future analysis. 
  #This is for the fall/summer recordings when air temperature doesn't get much below 0°C
  #This needs to be run prior to filtering accelerometer data by date so that every 10 seconds are the same temperature.
  #Studd 2017
  m<-seq(min(axy$dtime), max(axy$dtime), by="10 sec")
  m<-data.frame(dtime=m, id=seq(1, length(m), 1))
  axy<- left_join(axy, m, by="dtime") 
  axy<-axy %>% mutate(id=zoo::na.locf(id))
  sum<-axy %>% group_by(id) %>% summarise(temp=mean(temp))
  sum<-sum %>% mutate(difftemp1=temp-lag(temp), difftemp2=lead(temp)-temp) %>% mutate(tempC=ifelse(((difftemp1>=1.3 | difftemp1<=-1.3) & (difftemp2>=1.3 | difftemp2<=-1.3)), (lead(temp)+lag(temp))/2, temp))
  sum<-select(sum, id, tempC)
  axy<-left_join(axy, sum, by="id")
  axy<-select(axy, -id)
  return(axy)
}

tempFilterWinter= function(axy){
  ## This function filters the temperature recordings for the spurlious events where a faulty temperature is recorded. It then replaces that temperature with the median temperature between the recording preceeding and following it.  I double checked that these faulty temperature never seem to occur in consequative recordings. Temperature sensor only seems to register a temperature every 10 seconds (probably power saving method), so every 10 recordings have the same temperature.  Once this is used, tempC variable is teh temperature that you will want to use for all future analysis. 
  #This is for the winter recordings when air temperature is consistently well below 0°C
  #This needs to be run prior to filtering accelerometer data by date so that every 10seconds are the same temperature.
  #Studd 2017
  m<-seq(min(axy$dtime), max(axy$dtime), by="10 sec")
  m<-data.frame(dtime=m, id=seq(1, length(m), 1))
  axy<- left_join(axy, m, by="dtime") 
  axy<-axy %>% mutate(id=zoo::na.locf(id))
  sum<-axy %>% group_by(id) %>% summarise(temp=mean(temp))
  sum<-sum %>% mutate(difftemp1=temp-lag(temp), difftemp2=lead(temp)-temp) %>% mutate(tempC=ifelse((difftemp1>=1.5 | difftemp1<=-5) & (difftemp2>=1.5 | difftemp2<=-5), (lead(temp)+lag(temp))/2, temp))
  sum[1,5]=sum[1,2]
  sum<-select(sum, id, tempC)
  axy<-left_join(axy, sum, by="id")
  axy<-select(axy, -id)
  return(axy)
}

axyprep= function(axy, SN){
  
  axy<- axy %>% mutate(Xs=runmed(X,k=91, "median"), Ys=runmed(Y,k=91, "median"), Zs=runmed(Z,k=91, "median"), temp11=runmed(tempC, 11))
  
  m<-seq(min(axy$dtime), max(axy$dtime), by="24 hour")-12*60*60
  m<-data.frame(time=m, id=seq(1, length(m), 1))
  axy$group<-m[match(axy$dtime, m$time), which(names(m)=="id")]
  axy$group[1]=1
  axy<-mutate(axy, group=zoo::na.locf(group)) 
  m <- axy %>%  group_by(group) %>% summarise(center=kmeansR(temp11))
  #m <- axy %>%  group_by(group) %>% summarise(center=histTemp(temp11))
  axy<-left_join(axy, m, by="group")
  
  axy<-axy %>% mutate(var=rollapply(data=temp, width=420, by=30, FUN=var,  align="center", fill=NA))
  axy$var[1]=0
  axy<-axy %>% mutate(var=zoo::na.locf(var))
  
  
  #m<-seq(min(axy$dtime), max(axy$dtime), by="2 min")
  #m<-data.frame(dtime=m, id=seq(1, length(m), 1))
  #axy<- left_join(axy, m, by="dtime") 
  #axy<-axy %>% mutate(id=zoo::na.locf(id))
  
  #  var<-axy%>% group_by(id) %>% summarise(var=var(temp))
  #  axy<-left_join(axy, var, by="id")
  #  axy<-select(axy, -group, -date, -id, -time)
  
  
  
  axy$Squirrel<-SN
  return(axy)
}

kmeansR<-function(a){
  b=2
  model=kmeans(a, b)
  r=(model$centers[1]+model$centers[2])/2
  return(r)
}

summaryStats= function (axy) {
  axy<-axy %>% mutate(Xd=X-Xs, Yd=Y-Ys, Zd=Z-Zs, odba=abs(Xd)+abs(Yd)+abs(Zd)) %>% mutate(Dx=rollapply(data=Xd, width=11, by=11, FUN=diffX,  align="left", fill=NA), Dy=rollapply(data=Yd, width=11, by=11, FUN=diffX,  align="left", fill=NA), Dz=rollapply(data=Zd, width=11, by=11, FUN=diffX,  align="left", fill=NA), Dodba=Dx+Dy+Dz)
  axy<-axy %>% mutate(Dodba=zoo::na.locf(Dodba))
  
  #maxYd at 4sec, Todba at 10sec
  
  axy<-axy %>% mutate(maxY=rollapply(data=abs(Yd), width=5, by=5, FUN=max,  align="left", fill=NA), Todba=rollapply(data=odba, width=11, by=11, FUN=sum,  align="left", fill=NA))
  axy<-axy %>% mutate(maxY=zoo::na.locf(maxY), Todba=zoo::na.locf(Todba))
  return (axy)
}

nestCorrection= function (axy){
  axy1<-axy %>% mutate(id=seq(1, nrow(axy), 1), ID=ifelse(lag(Nest)==Nest, NA, id)) 
  axy1$ID[1]=1
  axy1<- axy1 %>% mutate(ID=as.numeric(factor(zoo::na.locf(ID)))) %>% select(-id) 
  a<-axy1 %>% group_by(ID, Nest, Move) %>% summarise(num=n())
  a<-a %>% spread(Move, num)
  a$Moving[is.na(a$Moving)]=1
  a$notMoving[is.na(a$notMoving)]=1
  a$ratio=a$Moving/a$notMoving
  a<-a %>% mutate(Nest2=ifelse(Nest=="Nest" & ratio<=1, "Nest", "Out"))
  a<-a %>% select(ID, Nest, Nest2)
  axy1<-left_join(axy1, a, by=c("ID", "Nest"))
  return(axy1$Nest2)
}

diffX<-function(data){
  ab=0
  for (i in 2:11){
    ab[i]<-data[i]- data[i-1]
  }
  r=sum(abs(ab))
  return(r)
}

## Use this for the classification. 
behavClass2=function(axy){
  axy<-axy %>% mutate(Nest=ifelse(temp11>center, "Nest", "Out"), Move=ifelse(Dodba<=1.06, "notMoving", "Moving"), Feed=ifelse(Move=="notMoving", "notmoving", ifelse(Todba>=6.2, "Forage", "Feed")), Travel=ifelse(Feed!="Forage", Feed, ifelse(maxY<1.15, "Forage", "Travel")))
  ## if reclassifying the nest data with the nest correction. Run this.  If working on lactating females, you may want to skip this next step  
  axy<-axy %>% mutate(Nest2=nestCorrection(axy)) %>% mutate(All=ifelse(Nest2=="Nest" & Move=="Moving", "NestMove", ifelse(Nest2=="Nest" & Move=="notMoving", "NestNotMove", Travel)))
  return(axy)
}

####   Running on an axy file for example ----
axy<-fread("XX_Feb27_14.csv", col.names=c("date", "time", "X", "Y", "Z", "temp"))
axy$dtime<-dmy_hms(paste(axy$date, substr(axy$time, 1, 8), sep=" "))
axy$date<-dmy(axy$date)
axy<-tempFilterFall(axy)
# filter out dates that you don't need
#axy<-filter(axy, date>=ymd("2014-02-09"),  date<=ymd("2014-02-14"))

# run the data through the calculations of the summary stats needed for the classification
axy<-axyprep(axy, "J11")
axy<-summaryStats(axy)
#classify the data
axy<-behavClass2(axy)

#plot the data in an actogram
axy1$Time<-as.POSIXct(strptime(substr(axy1$dtime, 12, 19), "%H:%M:%S"))
p1<-ggplot(axy1, aes(Time, date))+geom_tile(aes(fill=factor(All)), data=axy1)                  
