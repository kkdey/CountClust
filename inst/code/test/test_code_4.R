K=4;
G=100;
N=200;

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

K <- 4
barplot(t(omega_true),col=2:(K+1),axisnames=F,space=0,border=NA,main="",las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
title(main=paste("Structure Plot of initial chosen topic proportions,k=",K))


theta_b_true <- rbind(rdirichlet(2, c(0.2, 0.2, rep(0.6/G,G-2))),
                           rdirichlet(2, c(rep(0.5/G,G-3),0.3,0.1,0.1)),
                           rdirichlet(2, c(rep(0.4/G, G/2 -2),0.2,0.1,rep(0.3/G, G/2))),
                           rdirichlet(2, c(rep(0.6/G,G-2), 0.2, 0.2)));

read_counts <- t(do.call(cbind,lapply(1:dim(omega_true)[1], function(x)
                                                        {
                                                            if(Label.Batch[x]==1)
                                                              out <- rmultinom(1,1000,prob=omega_true[x,]%*%theta_b_true[c(2,4,6,8),]);
                                                            if(Label.Batch[x]==2)
                                                              out <- rmultinom(1,1000,prob=omega_true[x,]%*%theta_b_true[c(1,3,5,7),]);
                                                            return(out)
                                                             }  )));

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


library(permute);
library("BioPhysConnectoR");
perm_set=rbind(1:K,allPerms(1:K));
diff=array(0,dim(perm_set)[1]);
for (p in 1:dim(perm_set)[1])
{
  temp=docweights_2[,perm_set[p,]];
  diff[p]=fnorm(temp,omega_true);
}

p_star=which(diff==min(diff));
docweights_2=docweights_2[,perm_set[p_star,]];

barplot(t(docweights_2),col=2:(K+1),axisnames=F,space=0,border=NA,main="",las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
title(main=paste("Structure Plot of initial chosen topic proportions,k=",K))


