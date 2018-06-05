
rm(list=ls())

library(reshape2)
library(plyr)

library(ggplot2)

setwd('U:\\Gottwald_Lab\\Shuo Zhang\\CA sampling efficacy')



W_dist=function(a1,a2) round(sqrt((a1[,1]-a2[1])^2+(a1[,2]-a2[2])^2),0)  #calc dist



########## 

acp=read.csv('Code\\R data\\All_ACP_samples_with_STR_ID.csv',stringsAsFactors = F)

plant=read.csv('Code\\R data\\ALL STR_plant samples.csv',stringsAsFactors = F)


hlb=read.csv('Data\\riskHLB.csv',stringsAsFactors = F)
risk=read.csv('Data\\STR_risk_and_sampling_number.csv',stringsAsFactors = F)


#### ACP processing

acp1=aggregate(acp$STR_ID,by=list(acp$Collected_,acp$Longitude,acp$Latitude,acp$STR_ID),length)
names(acp1)=c('Date','Long','Lat','STR_ID','Tree_Sam')

sel=which(acp1$STR_ID>0)

acp=acp1[sel,]

acp$Date=as.Date(acp$Date,format = '%m/%d/%Y')



write.csv(acp1,'ACP sample locations up to 01_2018.csv',row.names = F)




acp$Year=as.numeric(format(acp$Date,'%Y'))

STR_ACP=aggregate(acp$Tree_Sam,by=list(acp$STR_ID,acp$Year),length)

STR_ACP2=aggregate(acp$Tree_Sam,by=list(acp$STR_ID,acp$Year),sum)

STR_ACP=data.frame(STR_ACP,STR_ACP2$x)

names(STR_ACP)=c('STR_ID','Year','Loc_Sam','Tree_Sam')

library(reshape2)

STR_ACP_1=dcast(STR_ACP,STR_ID~Year,value.var='Loc_Sam')

names(STR_ACP_1)[2:7]=paste('LACP',12:17,sep='')


STR_ACP_2=dcast(STR_ACP,STR_ID~Year,value.var='Tree_Sam')

names(STR_ACP_2)[2:7]=paste('ACP',12:17,sep='')

summary(STR_ACP_1[,1]-STR_ACP_2[,1])

STR_ACP_Fin=data.frame(STR_ACP_1,STR_ACP_2[,2:7])


##### Plant sample processing



plant1=aggregate(plant$TOTAL,by=list(plant$DATE,plant$Lat_OWNER,plant$Long_OWNER,plant$PLSFILL_),sum)
names(plant1)=c('Date','Lat','Long','STR_ID','Tree_Sam')



sel=which(plant1$STR_ID>0)

plant=plant1[sel,]

plant$Date=as.Date(plant$Date,format = '%m/%d/%Y')


plant$Year=as.numeric(format(plant$Date,'%Y'))

STR_plant=aggregate(plant$Tree_Sam,by=list(plant$STR_ID,plant$Year),length)

STR_plant2=aggregate(plant$Tree_Sam,by=list(plant$STR_ID,plant$Year),sum)

STR_plant=data.frame(STR_plant,STR_plant2$x)

names(STR_plant)=c('STR_ID','Year','Loc_Sam','Tree_Sam')


STR_plant_1=dcast(STR_plant,STR_ID~Year,value.var='Loc_Sam')

names(STR_plant_1)[2:7]=paste('Ltree',12:17,sep='')


STR_plant_2=dcast(STR_plant,STR_ID~Year,value.var='Tree_Sam')

names(STR_plant_2)[2:7]=paste('tree',12:17,sep='')

summary(STR_plant_1[,1]-STR_plant_2[,1])

STR_plant_Fin=data.frame(STR_plant_1,STR_plant_2[,2:7])



STR_Fin=merge(STR_ACP_Fin,STR_plant_Fin,by.x='STR_ID',by.y='STR_ID',all=T)

STR_Fin[is.na(STR_Fin)]=0

temp=STR_Fin[,2:13]+STR_Fin[,14:25]

names(temp)=c(paste('Loc',12:17,sep=''),paste('Tree',12:17,sep=''))

STR_total=data.frame(STR_Fin$STR_ID,temp)
names(STR_total)[1]='STR_ID'


STR_total=STR_total[order(STR_total[,1]),]


write.csv(STR_total,'STR_sampling density 12_17.csv',row.names = F)



##### HLB

HLB=hlb[,c(2,11)]

HLB$Collected=as.Date(HLB$Collected,format = '%m/%d/%Y')


HLB$Year=as.numeric(format(HLB$Collected,'%Y'))

STR_HLB=aggregate(HLB$PLSFILL,by=list(HLB$PLSFILL,HLB$Year),length)

names(STR_HLB)=c('STR_ID','Year','HLB_cnt')

STR_HLB=dcast(STR_HLB,STR_ID~Year,value.var='HLB_cnt')


STR_HLB[is.na(STR_HLB)]=0

names(STR_HLB)[2:6]=paste('HLB',c(11,12,15,16,17),sep='')


write.csv(STR_HLB,'STR_HLB 12_17.csv',row.names = F)




#### Estimate Actual HLB incidence...

STR_all=read.csv('Code\\R data\\CA_STR_sampling_risk_HLB_all_info.csv',stringsAsFactors = F)

STR_all=STR_all[order(STR_all$PLSFILL_),]

HLB_prop=matrix(0,nrow(STR_all),ncol=5)

HLB17_risk=STR_all$RCT2016
HLB17_risk=HLB17_risk/1.5
HLB17_risk[HLB17_risk>1]=1
  
PCR_sensitivity=0.25


### part 1 for STRs without HLB detection

sel=which(STR_all$Loc15_17>0) 
    HLB_prop[sel,1]=1-exp(log(0.95)/(STR_all$Loc15_17[sel]*PCR_sensitivity))
    HLB_prop[sel,2]=1-exp(log(0.5)/(STR_all$Loc15_17[sel]*PCR_sensitivity))
    HLB_prop[sel,3]=1-exp(log(0.05)/(STR_all$Loc15_17[sel]*PCR_sensitivity))

    HLB_prop[,4]=STR_all$Loc15_17
    HLB_prop[,5]=HLB17_risk

### part 2 for STR with HLB detection updates
    
sel=which(STR_all$HLB_all>0)

temp=STR_all[sel,]
aaa=data.frame(temp$HLB15/temp$Loc15,temp$HLB16/temp$Loc16,temp$HLB17/temp$Loc17)

HLB_prop2=matrix(0,length(sel),ncol=4)

HLB_prop2[,4]=round(PCR_sensitivity/apply(aaa,1,max,na.rm=T))

ss=which(HLB_prop2[,4]<1)
HLB_prop2[ss,4]=1


T_prob=seq(0.000001,1,by=0.000001)

for(i in 1:nrow(HLB_prop2)) {
  
  temp=1-pbinom(1,HLB_prop2[i,4],prob=T_prob)
  HLB_prop2[i,1]=T_prob[min(which(temp>0.05))]
  HLB_prop2[i,2]=T_prob[min(which(temp>0.5))]
  HLB_prop2[i,3]=T_prob[min(which(temp>0.95))]
  print(i)
}

HLB_prop[sel,1:3]=HLB_prop2[,1:3]

HLB_prop[is.na(HLB_prop)]=1

HLB_prop=data.frame(STR_all$PLSFILL_,HLB_prop)

names(HLB_prop)=c('STR_ID','Min_PHLB','Mean_PHLB','Max_PHLB','LOc15_17','Risk_adj')



write.csv(HLB_prop,'STR_actual_HLB_prob.csv',row.names = F)



#### known HLB+ dispersal risk


Plant_ct=read.csv('Code\\R data\\Plant_samples_with_CT_less_than_40_01_2018.csv',stringsAsFactors = F)

Plant_ct$Collected_=as.Date(Plant_ct$Collected_,format='%m/%d/%Y')

Plant_ct$Year=as.numeric(format(Plant_ct$Collected_,'%Y'))


ggplot(data=Plant_ct[,],aes(Results,CtValue))+geom_boxplot(col='green')+geom_jitter(width = 0.25)+labs(x='',y='CT value')

Plant_ct$CT_index=5-as.numeric(cut(Plant_ct$CtValue,c(0,32,36,38.5,40)))




ACP_ct=read.csv('Code\\R data\\ACP_samples_with_CT_less_than_40_01_2018.csv',stringsAsFactors = F)

ACP_ct$Collected_=as.Date(ACP_ct$Collected_,format='%m/%d/%Y')

ACP_ct$Year=as.numeric(format(ACP_ct$Collected_,'%Y'))

sel=which(ACP_ct$Results !='Fail')

ggplot(data=ACP_ct[sel,],aes(Results,CtValue))+geom_boxplot(col='blue')+geom_jitter(width = 0.25)+labs(x='',y='CT value')

ACP_ct$CT_index=5-as.numeric(cut(ACP_ct$CtValue,c(0,32,36,38.5,40)))



HLB_ct=read.csv('Code\\R data\\HLB_finds_locations_up_to_01_2018.csv',stringsAsFactors = F)

HLB_ct$Collected=as.Date(HLB_ct$Collected,format='%m/%d/%Y')

HLB_ct$Year=as.numeric(format(HLB_ct$Collected,'%Y'))




STR=read.csv('Code\\R data\\STR_03_21_2018.csv',stringsAsFactors = F)



#### ACP ct risk
ACP_ctR=matrix(0,nrow=nrow(STR),ncol=6)

a1=as.matrix(STR[,10:11])

CT_w=c(0,1,1.5,2)

my_year=2012:2017

for (j in 1:length(my_year)) {
  sel1=which(ACP_ct$Year==my_year[j] & ACP_ct$Results !='POSITIVE')
  
  if(length(sel1)>0) {
    New_value=ACP_ct[sel1,]    #### only the ACP with LAS <=36
    
  
    for (i in 1:nrow(New_value)) {
      
      a2=as.numeric(New_value[i,27:28])
      
      temp=W_dist(a1,a2)/1000
      
      sel=which(temp<=16)
      if(length(sel)>0) {
        ACP_ctR[sel,j]=ACP_ctR[sel,j]+1/(1+(temp[sel])^(1.6))*CT_w[New_value$CT_index[i]]
      }
      
      if (i %% 100 ==0) {print(i) }
      
    }
    
  }
  
}

ACP_ctR=data.frame(ACP_ctR)

names(ACP_ctR)=paste('ACPct',12:17,sep='')



#### Tree ct risk
Tree_ctR=matrix(0,nrow=nrow(STR),ncol=6)

a1=as.matrix(STR[,10:11])

CT_w=c(0,1,1.5,2)

my_year=2012:2017

for (j in 1:length(my_year)) {
  sel1=which(Plant_ct$Year==my_year[j] & Plant_ct$Results !='POSITIVE')
  
  if(length(sel1)>0) {
    New_value=Plant_ct[sel1,]    #### only the ACP with LAS <=36
    
    
    for (i in 1:nrow(New_value)) {
      
      a2=as.numeric(New_value[i,24:25])
      
      temp=W_dist(a1,a2)/1000
      
      sel=which(temp<=16)
      if(length(sel)>0) {
        Tree_ctR[sel,j]=Tree_ctR[sel,j]+1/(1+(temp[sel])^(1.6))*CT_w[New_value$CT_index[i]]
      }
      
      if (i %% 100 ==0) {print(i) }
      
    }
    
  }
  
}

Tree_ctR=data.frame(Tree_ctR)

names(Tree_ctR)=paste('Treect',12:17,sep='')



#### Confirmed HLB+ location risk
HLB_ctR=matrix(0,nrow=nrow(STR),ncol=6)

a1=as.matrix(STR[,10:11])

CT_w=2

my_year=2012:2017

for (j in 1:length(my_year)) {
  sel1=which(HLB_ct$Year==my_year[j])
  
  if(length(sel1)>0) {
    New_value=HLB_ct[sel1,]    #### 
    
    
    for (i in 1:nrow(New_value)) {
      
      a2=as.numeric(New_value[i,68:69])
      
      temp=W_dist(a1,a2)/1000
      
      sel=which(temp<=16)
      if(length(sel)>0) {
        HLB_ctR[sel,j]=HLB_ctR[sel,j]+1/(1+(temp[sel])^(1.6))*CT_w
      }
      
      if (i %% 100 ==0) {print(i) }
      
    }
    
  }
  
}

HLB_ctR=data.frame(HLB_ctR)

names(HLB_ctR)=paste('HLBct',12:17,sep='')



Fin_ctR=data.frame(STR$PLSFILL_,HLB_ctR,Tree_ctR,ACP_ctR)


HLB_FRisk=rowSums(Fin_ctR[,-1])

Fin_ctR=data.frame(Fin_ctR,HLB_FRisk)



write.csv(Fin_ctR,'STR_HLB_CT_risk 2018.csv',row.names = F)



###### Healthy ACP locatoins


ACP_Loc=read.csv('Code\\R data\\ACP_sample_aggregate_per_location projected_01_2018.csv',stringsAsFactors = F)


ACP_Loc$Date=as.Date(ACP_Loc$Date,format = '%m/%d/%Y')


ACP_Loc$Year=as.numeric(format(ACP_Loc$Date,'%Y'))



ACP_Risk=matrix(0,nrow=nrow(STR),ncol=6)

a1=as.matrix(STR[,10:11])

my_year=2012:2017

for (j in 1:length(my_year)) {
  sel1=which(ACP_Loc$Year==my_year[j])
  
  if(length(sel1)>0) {
    New_value=ACP_Loc[sel1,]    #### 
    
    
    for (i in 1:nrow(New_value)) {
      
      a2=as.numeric(New_value[i,6:7])
      
      temp=W_dist(a1,a2)/1000
      
      sel=which(temp<=16)
      if(length(sel)>0) {
        ACP_Risk[sel,j]=ACP_Risk[sel,j]+1/(1+(temp[sel])^(1.6))*sqrt(ACP_Loc$Tree_Sam[i])
      }
      
      if (i %% 100 ==0) {print(i) }
      
    }
    
  }
  
}


ACP_Risk=data.frame(ACP_Risk)

names(ACP_Risk)=paste('ACP_ri',12:17,sep='')


ACP_Risk=data.frame(STR$PLSFILL_,ACP_Risk)


write.csv(ACP_Risk,'STR_Heathy_ACP_risk 2018.csv',row.names = F)



Total_R17=(1*STR_R2017$CT_R17+2.57*STR_R2017$ACP_R17+0.43*STR_R2017$PACP_R17+0.78*STR_R2017$Road_R17+3*STR_R2017$LAS_R17+0.25*STR_R2017$PackH_R17+0.25*STR_R2017$FM_R17+0.1*STR_R2017$MIR_R17+0.1*STR_R2017$Organ_R17)*STR_R2017$Tmin_R17*log(STR_R2017$RC_tree+1)



write.csv(aa,'Inconclusive CT plant.csv',row.names = F)



save.image('Code\\sample efficacy part 3.RData')





