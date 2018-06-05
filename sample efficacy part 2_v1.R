
rm(list=ls())

library(reshape2)
library(plyr)

setwd('U:\\Gottwald_Lab\\Shuo Zhang\\CA sampling efficacy')

########## combine plant data

dat=read.csv('Data\\All_ACP_samples_up_to20180108.csv',stringsAsFactors = F)
dat$Collected_Date=as.Date(dat$Collected_Date,format = '%m/%d/%Y')

dat$Year=format(dat$Collected_Date,'%Y')
dat$MMDD=substring(dat$Collected_Date,6)

dat.xy=unique(dat[,c(2:6)])
names(dat.xy)[1]='Date'

sel=which(nchar(dat$City)==1)

agg=aggregate(dat$PDRNumber,by=list(dat$Year,dat$MMDD,dat$Address,dat$City),length)
names(agg)[3:5]=c('Address','City','N')

agg$Date=as.Date(paste(agg$Group.1,agg$Group.2,sep = '-'))

agg.m=merge(agg[,3:6],dat.xy,by=c('Date','Address','City'),all.x=T)

write.csv(agg.m,'Data\\ACP_sample_up_to_20180108.csv',row.names = F)

acp=read.csv('Data\\riskACP.csv',stringsAsFactors = F)
plant=read.csv('Data\\Risk_plant.csv',stringsAsFactors = F)
hlb=read.csv('Data\\riskHLB.csv',stringsAsFactors = F)
risk=read.csv('Data\\STR_risk_and_sampling_number.csv',stringsAsFactors = F)

acp$Collected_Date=as.Date(acp$Collected_Date,format = '%m/%d/%Y')
acp$Year=format(acp$Collected_Date,'%Y')
acp$MMDD=substring(acp$Collected_Date,6)
acp$Month=format(acp$Collected_Date,'%m')

plant$Collected_Date=as.Date(plant$Collected_Date,format='%m/%d/%Y')
plant$MMDD=substring(plant$Collected_Date,6)

hlb.ID=unique(hlb$Address)
hlb$Collected=as.Date(hlb$Collected,format="%m/%d/%Y")
hlb$year=format(hlb$Collected,'%Y')
hlb$MMDD=substring(hlb$Collected,6)
hlb$Month=format(hlb$Collected,'%m')
  
hlb.info=hlb[!duplicated(hlb[,'Address']),c(4,6:8,11,15:21)]

str.hlb=unique(hlb.info$PLSFILL)

Table_4=hlb.info[,c(1,5:12)]

Table_4$Thoery=ceiling(Table_4$RC_tree*Table_4$Sam_rate/2)

agg.hlb=aggregate(hlb$Host,by=list(hlb$Address,hlb$year,hlb$MMDD),length) ## aggregate by address, year and date

agg.hlb=agg.hlb[order(agg.hlb[,1]),]
names(agg.hlb)=c('Address','Year','MMDD','N')

table.hlb=merge(agg.hlb,hlb.info[,1:5],by='Address')

table.hlb$Date=as.Date(paste(table.hlb$Year,table.hlb$MMDD,sep='-'))

######### Table 1-3 #########

Table_1=table.hlb

prev.year=3

Table_1$Prev_3year=rep(0,nrow(Table_1))
Table_1$Prev_2year=rep(0,nrow(Table_1))
Table_1$Prev_1year=rep(0,nrow(Table_1))

Table_2=Table_1

t.year=names(Table_1)[10:12]

Table_3=as.data.frame(matrix(NA,nrow(Table_1),15))
names(Table_3)[1:9]=names(Table_1)[1:9]

t3.acp= c('ACP_3','ACP_2','ACP_1')
t3.plant=c('Plant_3','Plant_2','Plant_1')
names(Table_3)[10:15]=c(t3.acp,t3.plant)
Table_3[,1:9]=Table_1[,1:9]

for (i in 1:nrow(Table_1)) {
  
  date.range=seq(Table_1[i,]$Date-365*prev.year,Table_1[i,]$Date,1)  #date range 3 year before HLB find
  
  date.int=c(date.range[365],date.range[365*2])
  
  sel.a=which(acp$PLSFILL==Table_1[i,]$PLSFILL & acp$Collected_Date %in% date.range)  ## select the one in the same STR and date range
  
  if(length(sel.a)>0) {
    
    temp.a=acp[sel.a,]
    
    temp.a$int.find=findInterval(temp.a$Collected_Date,date.int)
    
    agg.a=aggregate(temp.a$Address,by=list(temp.a$Address,temp.a$MMDD,temp.a$int.find),length)
    
    n.acp=aggregate(agg.a$Group.2,by=list(agg.a$Group.3),length)
    
    cl.a1=t.year[(n.acp$Group.1+1)]
    Table_1[i,cl.a1]=as.numeric(n.acp$x)
    
    acp_days=ddply(agg.a,~Group.3,summarise,n_days=length(unique(Group.2)))
    
    cl.a2=t3.acp[(n.acp$Group.1+1)]
    Table_3[i,cl.a2]=as.numeric(acp_days$n_days)
  }
  
  
  sel.p=which(plant$PLSFILL==Table_1[i,]$PLSFILL & plant$Collected_Date %in% date.range)
  
  if(length(sel.p)>0){
    
    
    temp.p=plant[sel.p,]
    
    temp.p$int.find=findInterval(temp.p$Collected_Date,date.int)
    
    agg.p=aggregate(temp.p$Address,by=list(temp.p$Address,temp.p$MMDD,temp.p$int.find),length)
    
    n.plant=aggregate(agg.p$Group.2,by=list(agg.p$Group.3),length)
    
    cl.p1=t.year[(n.plant$Group.1+1)]
    Table_2[i,cl.p1]=as.numeric(n.plant$x)
    
    plant_days=ddply(agg.p,~Group.3,summarise,n_days=length(unique(Group.2)))
    
    cl.p2=t3.plant[(n.plant$Group.1+1)]
    Table_3[i,cl.p2]=as.numeric(plant_days$n_days)
    
  }
  
  print(i)
  
}

Table.1=Table_1
Table.1$Prev_3year=rowSums(Table_1[,t.year])
Table.1$Prev_2year=rowSums(Table_1[,t.year[-1]])

Table.2=Table_2
Table.2$Prev_3year=rowSums(Table_2[,t.year])
Table.2$Prev_2year=rowSums(Table_2[,t.year[-1]])


######### Table 5 ##########

risk.th=10

Table_5=merge(Table.1,Table_4[,-c(3:4)],by=c('Address','PLSFILL'),all.x=T)

Table_5$Prev_3year=Table.1$Prev_3year+Table.2$Prev_3year
Table_5$Prev_2year=Table.1$Prev_2year+Table.2$Prev_2year
Table_5$Prev_1year=Table.1$Prev_1year+Table.2$Prev_1year

Table_5$Total_3year=Table_5$Prev_3year
Table_5$Total_2year=Table_5$Prev_2year
Table_5$Total_1year=Table_5$Prev_1year
Table_5$Year=as.numeric(Table_5$Year)

tt.year=names(Table_5)[19:21]
risk.year=names(Table_5)[13:17]

for (i in 1:nrow(Table_5)) {
  
  sel=which(Table_5[i,t.year]< Table_5[i,'Thoery'])
  
  if(length(sel)>0) {
    
    
    cl.name=paste('RCT',(Table_5[i,'Year']-(4-sel)),sep='')
    
    sel.cl=which(cl.name %in% risk.year)
     
    if(length(sel.cl)>0) {
      
      tt=ceiling(Table_5[i,'Thoery']*(1+pmin(1,Table_5[i,cl.name]/risk.th)))
      
      Table_5[i,tt.year[sel[sel.cl]]]=tt
      
    }
    
    
    
  }
  
}

save.image('Code\\sample efficacy part2_v1.RData')

write.csv(Table.1,'Output\\Table_1.csv',row.names = F)
write.csv(Table.2,'Output\\Table_2.csv',row.names = F)
write.csv(Table_3,'Output\\Table_3.csv',row.names = F)
write.csv(Table_4,'Output\\Table_4.csv',row.names = F)
write.csv(Table_5,'Output\\Table_5.csv',row.names = F)

######### plot ###########

library(ggplot2)
library(plotly)

############ HLB finds and AcP visits 


hlb.p=aggregate(hlb$Host,by=list(hlb$County,hlb$year),length)
names(hlb.p)=c('County','Year','N')

hlb.p$Year=as.factor(hlb.p$Year)

ggplot(hlb.p,aes(x=Year,y=N,fill=County))+geom_bar(position = 'dodge', stat = 'identity')+
  labs(x='Year',y='No. HLB+')+theme_bw()+theme(legend.position=c(.2,.8),text=element_text(size=16))




hlb_1yr_acp=hlb[,c(2,4,6,9:11,22)]    ## look back 1 year to see how many acp visits in each STR

hlb_1yr_acp$ACP_visits=rep(0,nrow(hlb_1yr_acp))

for (i in 1:nrow(hlb)) {
  
  date.range=seq(hlb[i,]$Collected-365,hlb[i,]$Collected,1)  #date range 1 year before HLB find
  
  sel=which(acp$PLSFILL==hlb[i,]$PLSFILL & acp$Collected_Date %in% date.range)
  
  if(length(sel)>0){
    
    temp=acp[sel,]
    
    agg1=aggregate(temp$Address,by=list(temp$Address,temp$MMDD,temp$Year),length)
    
    hlb_1yr_acp[i,]$ACP_visits=nrow(agg1)
    
  }
  
}


agg1=aggregate(hlb_1yr_acp$ACP_visits,by=list(hlb_1yr_acp$PLSFILL,hlb_1yr_acp$Year,hlb_1yr_acp$County),length)
agg2=aggregate(hlb_1yr_acp$ACP_visits,by=list(hlb_1yr_acp$PLSFILL,hlb_1yr_acp$Year,hlb_1yr_acp$County),mean)

n_acp=merge(agg1,agg2,by=c('Group.1','Group.2','Group.3'))
names(n_acp)=c('STR','Year','County','N_HLB','N_ACP')

n_acp$ACP_visits=round(n_acp$N_ACP/n_acp$N_HLB)

sel=which(n_acp$County %in% c('Los Angeles','Orange'))

ggplot(n_acp[sel,],aes(x=Year,y=ACP_visits,fill=County))+geom_boxplot()+labs(x='Year',y='No. Visits')+
  theme_bw()+theme(legend.position=c(.2,.8),text=element_text(size=16))

unique(acp$Year)

##### 2017 HLB+ in each month

hlb.p=aggregate(hlb$Host,by=list(hlb$County,hlb$year,hlb$Month),length)
names(hlb.p)=c('County','Year','Month','N')

sel=which(hlb.p$Year==2016 & hlb.p$County %in% c('Los Angeles','Orange'))

hlb.p$Month=as.factor(hlb.p$Month)

ggplot(hlb.p[sel,],aes(x=Month,y=N,fill=County))+geom_bar(position = 'dodge', stat = 'identity')+
  labs(x='Month',y='No. HLB+')+theme_bw()+theme(legend.position=c(.2,.8),text=element_text(size=16))


sel=which(hlb$year==2016 & hlb$County %in% c('Los Angeles','Orange'))

hlb_acp=hlb[sel,c(2,4,6,9:11,22,24)]    ## look back __ year to see how many acp visits in each STR

hlb_acp$ACP_visits=rep(0,nrow(hlb_acp))

for (i in 1:nrow(hlb_acp)) {
  
  date.range=seq(hlb_acp[i,]$Collected-365,hlb_acp[i,]$Collected,1)  #date range 1 year before HLB find
  
  sel=which(acp$PLSFILL==hlb_acp[i,]$PLSFILL & acp$Collected_Date %in% date.range)
  
  if(length(sel)>0){
    
    temp=acp[sel,]
    
    agg1=aggregate(temp$Address,by=list(temp$Address,temp$MMDD,temp$Year),length)
    
    hlb_acp[i,]$ACP_visits=nrow(agg1)
    
  }
  
}

head(hlb_acp)

agg1=aggregate(hlb_acp$ACP_visits,by=list(hlb_acp$PLSFILL,hlb_acp$year,hlb_acp$Month,hlb_acp$County),length)
agg2=aggregate(hlb_acp$ACP_visits,by=list(hlb_acp$PLSFILL,hlb_acp$year,hlb_acp$Month,hlb_acp$County),mean)

n_acp=merge(agg1,agg2,by=c('Group.1','Group.2','Group.3','Group.4'))
names(n_acp)=c('STR','Year','Month','County','N_HLB','N_ACP')
head(n_acp)

n_acp$ACP_visits=round(n_acp$N_ACP/n_acp$N_HLB)

n_acp1=aggregate(n_acp$ACP_visits,by=list(n_acp$Month,n_acp$County),mean)

names(n_acp1)=c('Month','County','ACP_visits')
head(n_acp1)

sel=which(hlb.p$Year==2016 & hlb.p$County %in% c('Los Angeles','Orange'))

dat.p=merge(n_acp1,hlb.p[sel,],by=c('County','Month'),all.x=T)

tmp=max(dat.p$N)/max(dat.p$ACP_visits/2)

###2017
ggplot(dat.p,aes(x=Month,y=N/tmp,fill=County))+geom_bar(position = 'dodge', stat = 'identity')+
  labs(x='Month',y='')+theme_bw()+theme(legend.position=c(.2,.8),text=element_text(size=16))+
  geom_line(aes(x=as.numeric(Month),y=ACP_visits,group=County,colour=County))+
  geom_point(aes(x=as.numeric(Month),y=ACP_visits,group=County,colour=County,shape=County),size=3,fill='white')+
  scale_y_continuous(name = 'ACP visits',sec.axis = sec_axis(trans = ~.*tmp,name='No. HLB+'))

####2016
ggplot(dat.p,aes(x=Month,y=N/tmp))+geom_bar(fill='#A4A4A4',stat = 'identity')+
  labs(x='Month',y='')+theme_bw()+
  geom_line(aes(x=Month,y=ACP_visits,group=County,colour=County))+
  geom_point(aes(x=Month,y=ACP_visits,colour=County),size=3)+
  scale_y_continuous(name = 'ACP visits',sec.axis = sec_axis(trans = ~.*tmp,name='No. HLB+'))+
  theme(legend.position="none",text=element_text(size=16))

######## ACP by county

t1=Table.1[,c(5,8:12)]

dat.p1=melt(t1,id.vars = c('County','PLSFILL','Date'))
dat.p1$County=as.factor(dat.p1$County)

dat.p1$log.N=log(dat.p1$value)
sel=which(dat.p1$log.N==-Inf)
dat.p1[sel,]$log.N=0

ggplot(dat.p1,aes(x=variable,y=value,fill=County))+geom_boxplot()+labs(x='Year',y='No. Visits')+
  theme_bw()+theme(text=element_text(size=16))

ggplot(dat.p1,aes(x=variable,y=log.N,fill=County))+geom_boxplot()+labs(x='Year',y='log(No. Visits)')+
  theme_bw()+theme(text=element_text(size=16))

ggplot(dat.p1,aes(x=variable,y=sqrt(value),fill=County))+geom_boxplot()+labs(x='Year',y='sqrt(No. Visits)')+
  theme_bw()+theme(text=element_text(size=16))

 
dat.p1.mean=ddply(dat.p1,c('County','variable'),summarise,
                  N=length(value),
                  mean=mean(value),
                  sd=sd(value),
                  se=sd/sqrt(N)
                  )

ggplot(dat.p1.mean,aes(x=variable,y=mean,fill=County))+geom_bar(position = 'dodge', stat = 'identity')+labs(x='Year',y='No. Visits')+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2,position = position_dodge(0.9))+
  theme_bw()+theme(text=element_text(size=16))

####### ACP vs Total

t2=Table.1[,-c(4:9)]

#sel=which(t2$Prev_1year<Table_5$Thoery)
#t2[sel,]$Prev_1year=ceiling((1+pmin(1,Table_5[sel,]$RCT2016/risk.th))*Table_5[sel,]$Thoery)

t2.long=melt(t2,id.vars = c('Address','Year','MMDD'))

#t2.long$variable=substring(t2.long$variable,6)

t3=Table_5[,c(1,3:4,19:21)]
t3.long=melt(t3, id.vars = c('Address','Year','MMDD'))

t3.long$variable=gsub('Total','Prev',t3.long$variable)

t23=merge(t2.long,t3.long,by=c('Address','Year','MMDD','variable'))

names(t23)[4:6]=c('Prev_year','ACP','Total')

dat.p2=melt(t23,id.vars = c('Address','Year','MMDD','Prev_year'))
dat.p2$Prev_year=as.factor(dat.p2$Prev_year)
dat.p2$Prev_year=factor(dat.p2$Prev_year,levels(dat.p2$Prev_year)[3:1])

dat.p2$log.N=log(dat.p2$value)
sel=which(dat.p2$log.N==-Inf)
dat.p2[sel,]$log.N=0

ggplot(dat.p2,aes(x=Prev_year,y=value,fill=variable))+geom_boxplot()+labs(x='Year',y='No. Visits')+
  theme_bw()+theme(text=element_text(size=16))

ggplot(dat.p2,aes(x=Prev_year,y=sqrt(value),fill=variable))+geom_boxplot()+labs(x='Year',y='sqrt(No. Visits)')+
  theme_bw()+theme(text=element_text(size=16))

ggplot(dat.p2,aes(x=Prev_year,y=log.N,fill=variable))+geom_boxplot()+labs(x='Year',y='log(No. Visits)')+
  theme_bw()+theme(text=element_text(size=16))

acp.n=colMeans(t2[,4:6])
total.n=colMeans(t3[,4:6])

t23.mean=data.frame(ACP=acp.n,Total=total.n,Year=t.year)

t23.mean=ddply(dat.p2,c('Prev_year','variable'),summarise,
               N=length(value),
               mean=mean(value),
               sd=sd(value),
               se=sd/sqrt(N)
               )

ggplot(t23.mean,aes(x=Prev_year,y=mean,fill=variable))+geom_bar(position = 'dodge', stat = 'identity')+labs(y='No. Visits')+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2,position = position_dodge(0.9))+
  theme_bw()+theme(text=element_text(size=16))


########## Risk and 1HLB/# visits ##########

dat.r=merge(Table.1[,-c(4:7,9)],Table_4[,-c(3:4)],by=c('Address','PLSFILL'))

sel=which(rowSums(dat.r[,5:7])==0)   ### select the place with on data

dat.r=dat.r[-sel,]

Risk=rep(NA,nrow(dat.r))

r.year=as.numeric(dat.r$Year)-1  ## change number here 1,2,3

for (i in 1:nrow(dat.r)) {
  
  tmp=paste('RCT',r.year[i],sep='')  
  
  if( tmp %in% risk.year ) {Risk[i]=dat.r[i,tmp]}
  
}


Risk_1=Risk
Risk_2=Risk
Risk_3=Risk

Risk.1=as.character(cut(Risk_1,breaks = c(-0.1,0.1,0.5,1,5,500),labels = c('Low','Mid-low','Mid','Mid-high','High')))
Risk.2=as.character(cut(Risk_2,breaks = c(-0.1,0.1,0.5,1,5,500),labels = c('Low','Mid-low','Mid','Mid-high','High')))
Risk.3=as.character(cut(Risk_3,breaks = c(-0.1,0.1,0.5,1,5,500),labels = c('Low','Mid-low','Mid','Mid-high','High')))



prop_1=ceiling(1/dat.r$Prev_1year*100)
sel=which(dat.r$Prev_1year<dat.r$Thoery)
prop_1[sel]=ceiling(1/((1+pmin(1,Risk_1[sel]/risk.th))*dat.r[sel,]$Thoery)*100)

summary(prop_1)
prop_2=ceiling(1/dat.r$Prev_2year*100)
summary(prop_2)
prop_3=ceiling(1/dat.r$Prev_3year*100)
summary(prop_3)

Y=cbind(Find=prop_1,NotFind=100-prop_1)

dat=data.frame(Y, Risk=Risk.1, risk=Risk_1)

model.log=glm(cbind(Find,NotFind)~Risk,family = binomial(),data = dat)

summary(model.log)

model.pois=glm(Prev_1year~Risk.1,family = 'poisson',data=dat.r)

library(VGAM)

model.vg=vglm(dat.r$Prev_1year~Risk.1,family = poissonff)

model.vg1=vglm(prop_1~Risk.1,tobit(Upper = 100))


dat$Proportion=dat$Find/(dat$Find+dat$NotFind)
dat$N_visits=dat.r$Prev_1year

sel=which(dat.r$Prev_1year<dat.r$Thoery)
dat[sel,]$N_visits=(1+pmin(1,Risk_1[sel]/risk.th))*dat.r[sel,]$Thoery

ggplot(dat,aes(x=risk,y=round(1/Proportion)))+geom_point() +
  theme_bw()+theme(text=element_text(size=16))+stat_smooth(method = 'glm',col='red')

ggplot(dat,aes(x=risk,y=N_visits))+geom_point() +
  theme_bw()+theme(text=element_text(size=16))+stat_smooth(method = 'glm',col='red')

ggplot(dat,aes(x=risk,y=Find))+geom_point()+
  stat_smooth(method = 'glm',col='red')


ggplot(dat,aes(x=Risk,y=N_visits))+geom_boxplot()+
  theme_bw()+theme(text=element_text(size=16))





