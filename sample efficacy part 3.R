
rm(list=ls())

library(reshape2)
library(plyr)

setwd('U:\\Gottwald_Lab\\Shuo Zhang\\CA sampling efficacy')

########## 

acp=read.csv('Data\\STR_ACP.csv',stringsAsFactors = F)
plant=read.csv('Data\\STR_plant.csv',stringsAsFactors = F)
hlb=read.csv('Data\\riskHLB.csv',stringsAsFactors = F)
risk=read.csv('Data\\STR_risk_and_sampling_number.csv',stringsAsFactors = F)

acp$Date=as.Date(acp$Date,format = '%m/%d/%Y')
acp$Year=format(acp$Date,'%Y')
acp$MMDD=substring(acp$Date,6)
acp$Month=format(acp$Date,'%m')

plant$Date=as.Date(plant$DATE,format='%m/%d/%Y')
plant$Year=format(plant$Date,'%Y')
plant$MMDD=substring(plant$Date,6)

hlb.ID=unique(hlb$Address)
hlb$Collected=as.Date(hlb$Collected,format="%m/%d/%Y")
hlb$year=format(hlb$Collected,'%Y')
hlb$MMDD=substring(hlb$Collected,6)
hlb$Month=format(hlb$Collected,'%m')
  
hlb.info=hlb[!duplicated(hlb[,'Address']),c(4,6:8,11,15:21)]

strID=unique(acp$PLSFILL_,plant$PLSFILL_)
year=sort(unique(c(acp$Year,plant$Year)))

result=matrix(NA,ncol = 4)
colnames(result)=c('ID','Year','ACP_Sample','Total_Sample')


for (i in 1:length(strID)) {
  
  for (j in 1:length(year)) {
    
    
    sel.a=which(acp$PLSFILL_==strID[i] & acp$Year==year[j])
    sel.p=which(plant$PLSFILL_==strID[i] & plant$Year==year[j])
    
    temp=c(strID[i],year[j],length(sel.a),(length(sel.a)+length(sel.p)))
    
    result=rbind(result,temp)
  }
  
  print(i)
}


save.image('Code\\sample efficacy part 3.RData')





