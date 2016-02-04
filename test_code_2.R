
## test code for batch correction

K=4;
G=100;
N=200;
alpha_true=matrix(rnorm((K)*G,0.5,1),nrow=(K)); ### the matrix of fixed effects
Label.Batch=c(rep(1,N/2),rep(2,N/2));
B=max(Label.Batch);
sigmab_true=2;
beta_true=matrix(0,B,G);       ###  the matrix of the random effect
for(g in 1:G)
{
  beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}
library(gtools)
T=10;
omega_true=matrix(rbind(rdirichlet(T*4,c(3,4,2,6)),rdirichlet(T*4,c(1,4,6,3)),
                        rdirichlet(T*4,c(4,1,2,2)),rdirichlet(T*4,c(2,6,3,2)),
                        rdirichlet(T*4,c(3,3,5,4))), nrow=N);
over_dis=0.5;
noise_true=matrix(0,N,G);
for(n in 1:N)
{
  noise_true[n,]=over_dis*rnorm(G,0,1);
}

###  generating the table
read_counts=matrix(0,N,G);
for(n in 1:N)
{
  for(g in 1:G)
  {
    mean=exp(omega_true[n,]%*%alpha_true[,g] +beta_true[Label.Batch[n],g]+noise_true[n,g]);
    read_counts[n,g]=rpois(1,mean);
  }
}

library(CountClust)

library(maptpx)

K <- 4
Topic_clus <- topics(read_counts, K, tol=0.01);
docweights <- Topic_clus$omega;

barplot(t(docweights),col=2:(K+1),axisnames=F,space=0,border=NA,main="",las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
title(main=paste("Structure Plot of initial chosen topic proportions,k=",K))

### We do batch correction ####

batch_counts <- BatchCorrectedCounts(read_counts, batch_lab = Label.Batch);

Topic_clus_2 <- topics(batch_counts, K, tol=0.01);
docweights_2 <- Topic_clus_2$omega;

barplot(t(docweights_2),col=2:(K+1),axisnames=F,space=0,border=NA,main="",las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
title(main=paste("Structure Plot of initial chosen topic proportions,k=",K))
