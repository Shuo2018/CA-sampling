
rm(list=ls())

library(reshape2)

setwd('U:\\Gottwald_Lab\\Shuo Zhang\\CA sampling efficacy')

acp=read.csv('Data\\riskACP.csv',stringsAsFactors = F)
plant=read.csv('Data\\Risk_plant.csv',stringsAsFactors = F)
hlb=read.csv('Data\\riskHLB.csv',stringsAsFactors = F)
risk=read.csv('Data\\STR_risk_and_sampling_number.csv',stringsAsFactors = F)

acp$Collected_Date=as.Date(acp$Collected_Date,format = '%m/%d/%Y')
acp$Year=format(acp$Collected_Date,'%Y')
acp$MMDD=substring(acp$Collected_Date,6)

plant$Collected_Date=as.Date(plant$Collected_Date,format='%m/%d/%Y')
plant$MMDD=substring(plant$Collected_Date,6)

risk$RCT2017=rep(NA,nrow(risk))
risk$RCT2018=rep(NA,nrow(risk))
risk$Theory=round(risk$RC_tree*risk$Sam_rate/2)

risk1=risk[,c(3,9:16)]

fin.risk=melt(risk1,id.vars = c('PLSFILL','Theory'))
names(fin.risk)=c('STR','Theory','Year','Risk')

fin.risk$ACP=rep(0,nrow(fin.risk))
fin.risk$Plant=rep(0,nrow(fin.risk))

######### Table 1 ##########

year=2012:2018

risk.th=10

str.id=unique(risk$PLSFILL)


for (i in 1:length(str.id)) {
  
  sel.p=which(plant$PLSFILL==str.id[i])   ## select plant surveyed in each STR_ID
  sel.a=which(acp$PLSFILL==str.id[i])     ## select acp surveyed in each STR_ID
  
  if(length(sel.p)>0){
    
    temp.p=plant[sel.p,]     
    
    agg.p=aggregate(temp.p$Address,by=list(temp.p$Year,temp.p$MMDD,temp.p$Address),length)  ## aggregate by address, year and date
    agg.p1=aggregate(agg.p$Group.2,by=list(agg.p$Group.1),length)  ##how many places surveyed each year. (same place surveyed on 2 differect days considered as 2)
    
    for (k in 1:nrow(agg.p1)) {
      
      year.p=paste('RCT',agg.p1[k,1],sep='')
      sel.p.k=which(fin.risk$STR==str.id[i] & fin.risk$Year==year.p)
      
      if(length(sel.p.k)>0) {fin.risk[sel.p.k,]$Plant=agg.p1[k,2]}
      
    }
    
  }
  
  if(length(sel.a)>0){
    
    temp.a=acp[sel.a,]
    agg.a=aggregate(temp.a$Address,by=list(temp.a$Year,temp.a$MMDD,temp.a$Address),length)  ## aggregate by address, year and date
    agg.a1=aggregate(agg.a$Group.2,by=list(agg.a$Group.1),length)
    
    for (j in 1:nrow(agg.a1)) {
      
      year.a=paste('RCT',agg.a1[j,1],sep='')
      sel.a.j=which(fin.risk$STR==str.id[i] & fin.risk$Year==year.a)
      
     if(length(sel.a.j)>0) { fin.risk[sel.a.j,]$ACP=agg.a1[j,2]}
     
    }
    
  }
  
  print(i)
  
}

fin.risk$Total=(fin.risk$ACP+fin.risk$Plant)*pmin(1,fin.risk$Risk/risk.th)

write.csv(fin.risk,'Output\\table_1.csv',row.names = F)


######### Table 2 ACP only #######

library(geosphere)

acp2=read.csv('Data\\All_ACP_samples_up_to20180108.csv',stringsAsFactors = F)

acp2$Collected_Date=as.Date(acp2$Collected_Date,format = '%m/%d/%Y')
acp2$Year=format(acp2$Collected_Date,'%Y')
acp2$MMDD=substring(acp2$Collected_Date,6)

hlb.ID=unique(hlb$Address)
hlb$Collected=as.Date(hlb$Collected,format="%m/%d/%Y")
hlb$year=format(hlb$Collected,'%Y')
hlb$MMDD=substring(hlb$Collected,6)

agg.hlb=aggregate(hlb$Host,by=list(hlb$Address,hlb$year,hlb$MMDD),length) ## aggregate by address, year and date
agg.hlb1=aggregate(agg.hlb$Group.3,by=list(agg.hlb$Group.1,agg.hlb$Group.2),length) ## how many places surveyed each year. (same place surveyed on 2 differect days considered as 2)

agg.wide=dcast(agg.hlb1,Group.1~Group.2,value.var = 'x')

result=matrix(0,length(hlb.ID),(length(year)+1))

colnames(result)=c('Address',year)
result[,1]=hlb.ID

for(i in 1:length(hlb.ID)){
  
  hlb.xy=hlb[which(hlb$Address==hlb.ID[i])[1],c('Lon','Lat')]
  
  tmp.dist=distGeo(acp2[,c('Longitude','Latitude')],hlb.xy)
  
  sel=which(tmp.dist<=1609)  ## select distance within 1 mile
  
  if(length(sel)>0) {
    
    tmp.acp=acp2[sel,]
  
    agg.tmp=aggregate(tmp.acp$Address,by=list(tmp.acp$Year,tmp.acp$MMDD,tmp.acp$Address),length) ## aggregate by year, date and address
    agg.tmp1=aggregate(agg.tmp$Group.2,by=list(agg.tmp$Group.1),length)   ## how many places surveyed each year. (same place surveyed on 2 differect days considered as 2)
    
    sel.r=which(colnames(result) %in% agg.tmp$Group.1 )
    sel.a=which(agg.tmp1$Group.1 %in% colnames(result))
    
    result[i,sel.r]=agg.tmp1[sel.a,2]
    
  }
  
  print(i)
  
}

names(agg.wide)[1]='Address'

temp=as.matrix(agg.wide)
temp[which(is.na(temp)=='TRUE')]=0

temp=as.data.frame(temp)

result.m=merge(temp,result,by = 'Address')

hlb.ll=hlb[!duplicated(hlb[,'Address']),c('Address','Lat','Lon')]

result.m=merge(result.m,hlb.ll,by='Address')

write.csv(result.m,'Output\\table_2_xy.csv',row.names = F)

save.image('Code\\sample efficacy part 2.RData')

##########

library(ggplot2)

table2=read.csv('Output\\table_2.csv',stringsAsFactors =F )

ff=rep(0,nrow(result.m))

yr=c(2011:2012,2015:2017)

for (i in 1:nrow(table2)) {
  
  sel=which(table2[i,2:6]>0)
  
  ff[i]=yr[sel[1]]
  
}

table2$First.find=ff

dat.p=table2[,c(1,7:13,16)]

dat.p[,1]=1:nrow(dat.p)
names(dat.p)[1]='ID'

dat.p1=melt(dat.p,id.vars = c('ID','First.find'))

dat.p1$variable=substring(dat.p1$variable,5)
dat.p1$First.find=as.factor(dat.p1$First.find)

ggplot(dat.p1,aes(x=variable,y=value,fill=First.find))+geom_boxplot()+labs(x='Year',y='Samples')




