
setwd('U:\\Gottwald_Lab\\Shuo Zhang\\HLB sampling efficacy')


sample.d=c(seq(1,9,1),seq(10,100,5))/100
sample.n=c(seq(1,9,1),seq(10,100,5))





Sam=c(100,250,500,1000,5000)
D_prob=c(0.01,0.05,seq(0.1,0.5,0.1))

S_change=c(seq(1,9,1),seq(10,100,5))/100


result=matrix(0,nrow=length(S_change),ncol=length(D_prob)*length(Sam))


for (k in 1:length(D_prob)) {
  for (j in 1:length(Sam)) {
    for (i in 1:length(S_change)) {
      result[i,(k-1)*length(Sam)+j]=1-pbinom(0,round(Sam[j]*S_change[i]),D_prob[k])
    }
  }
  
}


colnames(result)=rep(Sam,length(D_prob))
rownames(result)=S_change

Sam=c(seq(1,5,1),seq(10,100,10),250,500)
D_prob=c(0.0001,0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.25,0.5)



result=matrix(0,nrow=length(Sam),ncol=length(D_prob))


for (k in 1:length(D_prob)) {
  for (j in 1:length(Sam)) {
   
      result[j,k]=1-pbinom(0,Sam[j],D_prob[k])
    }
}
  
colnames(result)=D_prob
rownames(result)=Sam



write.csv(result,file="HLB Sampling efficacy_2.csv")


