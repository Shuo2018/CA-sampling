
rm(list=ls())

library(ggplot2)
library(grid)
library(ggrepel)
library(RColorBrewer)
library(reshape2)

load("U:/Gottwald_Lab/Shuo Zhang/CA sampling efficacy/Code/HLB cluster analysis.RData")

## https://gis.stackexchange.com/questions/17638/how-to-cluster-spatial-data-in-r
## https://gis.stackexchange.com/questions/64392/find-clusters-of-points-based-distance-rule

library(sp)
library(rgdal)
library(geosphere)

############## visualize all the spatial points ##########

setwd('U:\\Gottwald_Lab\\Shuo Zhang\\CA sampling efficacy\\Data')

mydata=read.csv('HLB_finds_01_2018.csv',stringsAsFactors = F)

names(mydata)

dat=mydata[,-c(3,7:8,13:15,17:22)]

head(dat)
str(dat)

dat$Collected=as.Date(dat$Collected,'%m/%d/%Y')

dat$Year=format(dat$Collected,'%Y')
dat$Month=format(dat$Collected,'%m')

dat=dat[order(dat$FID_HLB_fi),]

dat.xy=dat[,c(11:12)]
row.names(dat.xy)=dat$FID_HLB_fi

dist.fid=as.matrix(round(dist(dat.xy)))

sort.dist=apply(dist.fid,1,FUN=sort)

min.dist=sort.dist[2,]

sel1=which(min.dist>800)
sel2=which(dat$Year %in% c(2011,2012))


dat.p=dat[,c(1,6,11:14)]

dat.p$iso=rep(0,nrow(dat.p))
dat.p[sel1,]$iso=1

cl=brewer.pal(9,'Set1')


x11()
ggplot(dat.p,aes(x=X_ref,y=Y_ref,color=as.factor(Year),shape=HostType))+
  geom_point(aes(size=as.factor(iso)))+
  scale_color_manual(values = cl[1:5])+
  scale_shape_manual(values = c(3,2))+
  scale_size_manual(values = c(1.5,3))+
  theme_bw()+
  #theme(legend.position="none")
  theme(legend.position=c(.8,.4))


########## cluster spatial points ##############

dat1=dat[-c(sel1,sel2),]
dat.xy1=dat1[,11:12]
row.names(dat.xy1)=dat1[,]$FID_HLB_fi

unique(dat1$Year)

dat.p=dat1[,c(1,6,11:14)]


cl=brewer.pal(5,'Set1')

sel=which(dat1$Year==2017)

dat.p=dat1[sel,]

ggplot(dat.p,aes(x=X_ref,y=Y_ref,color=HostType,shape=HostType))+
  geom_point(size=3)+
  #xlim(174000,250000)+
  #ylim(-480000,-430000)+
  #geom_text(aes(x=X_ref,y=Y_ref+c(1:24)*10,label=c(1:24)))+   # FID_HLB_fi[1]
  #annotate("text",x=x_a,y=y_a,label=c(1:24),size=3)+
  scale_color_manual(values = cl[1:2])+
  scale_shape_manual(values = c(3,2))+
  scale_size_manual(values = 2)+
  theme_bw()+
  #theme(legend.position="none",text=element_text(size=16))
  theme(legend.position=c(.7,.8),text=element_text(size=16))



sel=which(dat1$Year==2015)

xy=dat.xy1[sel,]

mdist=as.matrix(round(dist(xy)))


hc=hclust(as.dist(mdist),method = 'single')

d=seq(800,100,-100)

for (i in d) {
  
  temp=cutree(hc,h=i)
  cl.name=paste('clust_',i,sep = '')
  
  xy=data.frame(xy,temp)
  names(xy)[ncol(xy)]=cl.name
  
}

#xy

xy$HostType=dat1[sel,]$HostType
xy.2017p=xy

dat.p=xy
dat.p=xy_1
unique(dat.p$clust_800)

ggplot(dat.p,aes(x=X_ref,y=Y_ref,color=as.factor(clust_800)))+
  geom_point(size=3)+
  #scale_color_manual(values = cl[1:5])+
  theme_bw()+
  theme(text=element_text(size=16))    #######  legend.position=c(.7,.4),
  theme(legend.position="none",text=element_text(size=16))
  


tb=table(xy$clust_800)
tb
cl.n=which(as.numeric(tb)>20)
cl.n

sel=which(xy$clust_800==1)
sel=which(xy$clust_800==20 & xy$HostType=='Plant')

xy_1=xy[sel,c(1:2,11)]   ### go back to line 99-117 set xy=xy_1 then run line 144-155

mdist=as.matrix(round(dist(xy_1[,1:2])))

hc=hclust(as.dist(mdist),method = 'complete')

d=seq(800,100,-100)

for (i in d) {
  
  temp=cutree(hc,h=i)
  cl.name=paste('clust_',i,sep = '')
  
  xy_1=data.frame(xy_1,temp)
  names(xy_1)[ncol(xy_1)]=cl.name
  
  
  
}


dat.p=xy_1

ggplot(dat.p,aes(x=X_ref,y=Y_ref,color=as.factor(clust_300),shape=HostType))+
  geom_point(size=3)+
  #scale_color_manual(values = cl[1:5])+
  scale_shape_manual(values = c(3,2))+
  theme_bw()+
  theme(legend.position=c(.6,.7),text=element_text(size=16))    #######  legend.position=c(.7,.4),
  theme(legend.position="none",text=element_text(size=16))

temp_xy=data.frame(ID=row.names(xy),xy[,3:10])
temp_xy=data.frame(ID=row.names(xy_1),xy_1[,4:11])

dat.r=melt(temp_xy,id.vars = 'ID')

sel=which(dat.r$value==1)

dat.r$Find=rep(NA,nrow(dat.r))
dat.r[sel,]$Find='Yes'
dat.r[-sel,]$Find='No'

dat.r$value=as.factor(dat.r$value)

ggplot(dat.r,aes(x=variable,fill=Find))+geom_bar(position = 'fill')+
  theme_bw()+theme(text=element_text(size=16))





########## geom_path ######


### geom_path

tb=table(xy.all$clust_800)
tb


s.id=which(as.numeric(tb)>1)

j=1

sel=which(xy.all$clust_800==s.id[j])


dat.p=dat1[sel,]

dat.p$Days=as.numeric(dat.p$Collected-rep(dat.p[1,]$Collected,nrow(dat.p)))

sel=which(dat.p$Days==0)
p.size=rep(2,nrow(dat.p))
p.size[sel]=5
p.cl=rep(cl[1],nrow(dat.p))
p.cl[sel]=cl[2]

ggplot(dat.p,aes(x=X_ref,y=Y_ref))+
  geom_point(color=p.cl,size=p.size) +
  #geom_segment(aes(xend=c(tail(X_ref,n=-1),NA),yend=c(tail(Y_ref,n=-1),NA)),arrow = arrow(length = unit(0.1,'inches')),size=1,color=cl[2])+
  geom_text_repel(aes(X_ref,Y_ref,label=as.character(Collected)),size=5)+
  theme_classic(base_size = 16)


ggplot(dat.p,aes(x=X_ref,y=Y_ref,color=as.factor(Year)))+
  geom_point(size=3) +
  #geom_segment(aes(xend=c(tail(X_ref,n=-1),NA),yend=c(tail(Y_ref,n=-1),NA)),arrow = arrow(length = unit(0.1,'inches')),size=1,color=cl[2])+
  #geom_text_repel(aes(X_ref,Y_ref,label=FID_HLB_fi),size=5)+
  theme_classic(base_size = 16)+
  labs(title=paste('cluster',s.id[j]))

d=seq(100,800,50)

fin=matrix(0,length(s.id),length(d))
colnames(fin)=d
rownames(fin)=s.id

result1=NULL  ## clusters with more than 1 first points
result2=NULL  ## clusters with more than 6 months points
n_points=NULL ## number of points in each cluster

for (j in 1:length(s.id)) {
  
  sel=which(xy.all$clust_800==s.id[j])
  
  
  dat.p=dat1[sel,]
  
  dat.p$Days=as.numeric(dat.p$Collected-rep(dat.p[1,]$Collected,nrow(dat.p)))
  
  xy=unique(xy.all[sel,1:2])
  
  mdist=as.matrix(round(dist(xy)))
  
  n_points[j]=nrow(xy)
  
  if(max(dat.p$Days)<=180){
    
    sel.first=which(dat.p$Days==0)
    
    if(length(sel.first)>1 | nrow(xy)<2){
      
      result1=c(result1,s.id[j])
      
    }
    
    else{
      
      for (i in 1:length(d)) {
        
        temp1=which(mdist[,1]<d[i])
        temp2=1:nrow(mdist)
        temp3=temp1
        
        while (length(temp2)>length(temp1)) {
          
          temp1=temp3
          new.dist=as.matrix(mdist[temp3,])
          
          temp2=which(apply(new.dist, 2, min)<d[i])
          
          temp3=temp2
          
          print(paste('temp2=',length(temp2),', temp1=',length(temp1)))
          
        }
        
        fin[j,i]=round(length(temp3)/nrow(mdist)*100)
        print(d[i])
        
      }
      
      
      
      
      
    }
    
    
    
    
  }
  
  else{
    
    result2=c(result2,s.id[j])
  }
  
  print(paste('cluster',s.id[j],'done'))
}




for (i in 1:length(result2)) {
  
  sel=which(xy.all$clust_800==result2[i])
  
  
  dat.p=unique(dat1[sel,-c(1,6:10)])
  
  dat.p$Days=as.numeric(dat.p$Collected-rep(dat.p[1,]$Collected,nrow(dat.p)))
  
  dat.p$ct=c(rep(1,3),rep(2,3))   ## i=1: c(rep(1,6),rep(2,6),rep(3,2),rep(4,12),rep(5,1)) 
                                           ## i=2: c(rep(1,2),rep(2,2),rep(3,10)) 
                                           ## i=3: c(rep(4,1),rep(1,2),rep(2,2)) 
                                           ## i=4: c(rep(1,3),rep(2,3))
  ct=1:2
  
  fin_2=matrix(0,length(ct),length(d))
  
  for (j in 1:length(ct)) {
    
    sel.ct=which(dat.p$ct==ct[j])
    
    xy=dat.p[sel.ct,c('X_ref','Y_ref')]
    
    mdist=as.matrix(round(dist(xy)))
    
    for (k in 1:length(d)) {
      
      temp1=which(mdist[,1]<d[k])
      temp2=1:nrow(mdist)
      temp3=temp1
      
      while (length(temp2)>length(temp1)) {
        
        temp1=temp3
        new.dist=as.matrix(mdist[temp3,])
        
        temp2=which(apply(new.dist, 2, min)<d[k])
        
        temp3=temp2
        
        print(paste('temp2=',length(temp2),', temp1=',length(temp1)))
        
      }
      
      fin_2[j,k]=round(length(temp3)/nrow(mdist)*100)
    
    
  }
  
}

  
}

fin_2_1=fin_2
fin_2_2=fin_2
fin_2_3=fin_2
fin_2_4=fin_2

fin_2=rbind(fin_2_1,fin_2_2,fin_2_3,fin_2_4)

for (i in 1:length(result1)) {
  
  sel=which(xy.all$clust_800==result1[i])
  
  dat.p=unique(dat1[sel,-c(1,6:10)])
  
  dat.p$Days=as.numeric(dat.p$Collected-rep(dat.p[1,]$Collected,nrow(dat.p)))
  
  sel.first=which(dat.p$Days==0)
  
  temp=dat.p[-sel.first,]
  
  fin_1=matrix(0,length(sel.first),length(d))
  
  for (j in 1:length(sel.first)) {
    
    temp_new=rbind(dat.p[sel.first[j],],temp)
    
    xy=temp_new[,c('X_ref','Y_ref')]
    #xy=dat.p[,c('X_ref','Y_ref')]
    
    mdist=as.matrix(round(dist(xy)))
    
    for (k in 1:length(d)) {
      
      temp1=which(mdist[,1]<d[k])
      temp2=1:nrow(mdist)
      temp3=temp1
      
      while (length(temp2)>length(temp1)) {
        
        temp1=temp3
        new.dist=as.matrix(mdist[temp3,])
        
        temp2=which(apply(new.dist, 2, min)<d[k])
        
        temp3=temp2
        
        print(paste('temp2=',length(temp2),', temp1=',length(temp1)))
        
      }
      
      fin_1[j,k]=round(length(temp3)/nrow(mdist)*100)
      
      
    }
    
  }
  
}


fin_1_1=fin_1
fin_1_2=fin_1[2,]  ## cluster_11 2 first points are in the same location
fin_1_3=fin_1[1,]  ## cluster_12 2 first points are in the same location
fin_1_4=fin_1      ## cluster_17 3 first points have the same result
fin_1_5=fin_1[1,]  ## cluster_20 2 first points are in the same location
fin_1_6=fin_1[1,]  ## cluster_22 3 first points are in the same location

fin_1=rbind(fin_1_1,fin_1_2,fin_1_3,fin_1_4[1,],fin_1_5,fin_1_6)

fin[which(rownames(fin) %in% result1),]=fin_1
fin=fin[-which(rownames(fin) %in% result2),]
fin=rbind(fin,fin_2)

save.image('HLB cluster analysis.RData')

ID=1:nrow(fin)
df=data.frame(ID,fin)

df.m=melt(df,id.vars = 'ID')

df.m$variable=substring(df.m$variable,2)

p=ggplot(data = df.m,aes(x=variable,y=value))+geom_boxplot(colour=cl[2],fill=cl[5])+
  labs(x='Distance(m)',y='Percentage of coverage')+
  theme_bw()+theme(text=element_text(size=16))

new_df=aggregate(df.m$value,by=list(df.m$variable),mean)


p+geom_line(data=new_df,aes(x=as.factor(new_df$Group.1),y=new_df$x,group=1),size=1.5,col='red')




######## ACP HLB locations ##########

W_dist=function(a1,a2) round(sqrt((a1[,1]-a2[1])^2+(a1[,2]-a2[2])^2),0)  #calc dist

dat2=unique(dat[,c(2,6,11:12)])

sel=which(dat2$HostType=='ACP')

dat.acp=dat2[sel,]
dat.plant=dat2[-sel,]

m_dist=matrix(0,nrow(dat.plant),nrow(dat.acp))
m_days=matrix(0,nrow(dat.plant),nrow(dat.acp))

L=list()
N=0
for (i in 1:nrow(dat.acp)) {
  
  m_dist[,i]=W_dist(as.matrix(dat.plant[,3:4]),as.matrix(dat.acp[i,3:4]))
  m_days[,i]=as.numeric(dat.acp[i,]$Collected-dat.plant$Collected)
  
  sel=which(m_dist[,i]==0 & m_days[,i]==0)
  N[i]=length(sel)
  L[[i]]=sel
}


d=seq(50,1500,50)

result3=0

for (i in 1:length(d)) {
  
  temp=NA
  
  for (j in 1:ncol(m_dist)) {
  
    sel=which(m_dist[,j]<=d[i] & m_dist[,j]!=0)
    
    if(length(sel)>0){temp[j]=1} else{temp[j]=0}
    
    
  }
  
  result3[i]=round(sum(temp)/length(temp),2)*100
  
}

dat.p=data.frame(Distance=d,Percentage=result3)



sel_acp=0

for (i in 1:nrow(dat.acp)) {
  
  sel1=which(m_days[,i]<0)
  
  sel3=which(m_dist[,i]>0)
  temp=m_dist[sel3,i]
  
  sel2=which.min(temp)
  
  if(sel3[sel2] %in% sel1) {sel_acp[i]=0} else {sel_acp[i]=1}
}


sel=which((N+sel_acp)==1)


result4=0

m_dist1=m_dist[,sel]

for (i in 1:length(d)) {
  
  temp=NA
  
  for (j in 1:ncol(m_dist1)) {
    
    sel=which(m_dist1[,j]<=d[i] & m_dist1[,j]!=0)
    
    if(length(sel)>0){temp[j]=1} else{temp[j]=0}
    
    
  }
  
  result4[i]=round(sum(temp)/length(temp),2)*100
  
}

dat.p=data.frame(dat.p,Percentage=result4)
names(dat.p)[2:3]=c('All','Subset')

dat.m=melt(dat.p,id.vars = 'Distance')

dat.m$Distance=as.factor(dat.m$Distance)

ggplot(dat.m,aes(x=Distance,y=value,color=variable,group=variable))+
  geom_point()+geom_line()+
  labs(y='Percentage')+
  ylim(50,100)+
  theme_classic(base_size = 16)+
  theme_bw()


sel=which(sel_acp==1)

m_days1=m_days[,sel]

result5=0  ## how many days find the neareast tree
result6=0  ## how far away to the neareast tree

for (i in 1:ncol(m_days1)) {
  
  sel=which(m_dist1[,i]>0 & m_days1[,i]>0)
  
  temp=m_dist1[sel,i]
  
  sel1=which.min(temp)
  
  result6[i]=temp[sel1]
  result5[i]=m_days1[sel[sel1],i]
  
}

sel=which(sel_acp==1)

dat.p=data.frame(ID=sel,month=ceiling(result5/30))

ggplot(dat.p,aes(x=month))+geom_histogram(bins = 30,color='white',fill=cl[2])+
  theme_classic(base_size = 16)+
  theme_bw()






