# ABC model to estimate embryonic patch sizes in the human colon
# Tim Coorens - February 2021 
options(stringsAsFactors = F)

# Load in observed data
physical_dists=readRDS("physical_dists_colon.Rdata")
colon_dists=readRDS("genetic_dists_colon.Rdata")

colon_dists_vec=unlist(colon_dists)
colon_dists_vec=as.numeric(colon_dists_vec[!is.na(colon_dists_vec)])

physical_dists_vec=unlist(physical_dists)
physical_dists_vec=as.numeric(physical_dists_vec[!is.na(physical_dists_vec)])

freq_dist=table(physical_dists_vec)

# Run the simulation
sim_iter=50000
limit=as.numeric(names(freq_dist[freq_dist>10])[sum(freq_dist>10)])

prob_vec_obs_all=prob_vec_obs_all_total=matrix(0,nrow=sim_iter,ncol=limit)
n_iter=sum(physical_dists_vec%in%c(1:limit))
radius_vec=rep(NA,sim_iter)
x_max=max(physical_dists_vec)

freqs=table(physical_dists_vec)
prob_vec=rep(0,x_max)
names(prob_vec)=1:x_max
prob_vec[names(freqs)]=freqs

for(k in 1:sim_iter){
  radius=sample(seq(2,20,by=0.01),1)
  dist_vec=in_circle_vec=rep(NA,n_iter)
  
  for(n in 1:n_iter){
    crypt1_x=sample(seq(-radius,radius,by=1),1)
    y_lim=sqrt(radius^2-crypt1_x^2)
    crypt1_y=sample(seq(-y_lim,y_lim,by=1),1)
    
    dist_vec[n]=dist=sample(1:limit,1,prob=prob_vec[1:limit])
    angle=sample(seq(1,2*pi,by=0.01),1)
    crypt2_x=crypt1_x+cos(angle)*dist
    crypt2_y=crypt1_y+sin(angle)*dist
    
    
    in_circle_vec[n]=sqrt(crypt2_x^2+crypt2_y^2)<=radius
  }
  freqs_sharing_obs=table(dist_vec[in_circle_vec])
  freqs_all_obs=table(dist_vec)
  prob_vec_obs=prob_vec_obs_total=rep(0,limit)
  names(prob_vec_obs)=names(prob_vec_obs_total)=1:limit
  prob_vec_obs[names(freqs_sharing_obs)]=freqs_sharing_obs
  prob_vec_obs_total[names(freqs_all_obs)]=freqs_all_obs
  prob_vec_obs_all[k,]=prob_vec_obs
  prob_vec_obs_all_total[k,]=prob_vec_obs_total
  
  radius_vec[k]=radius
  if (k%%100==0){
    print(k)
  }
}

freqs_sharing2=table(physical_dists_vec[colon_dists_vec>15]) #Select crypts that share more than 15 SNVs
prob_vec_sharing2=rep(0,x_max)
names(prob_vec_sharing2)=1:x_max
prob_vec_sharing2[names(freqs_sharing2)]=freqs_sharing2

prob_vec_obs_all=as.data.frame(prob_vec_obs_all)
prob_vec_obs_all_total=as.data.frame(prob_vec_obs_all_total)
colnames(prob_vec_obs_all)=colnames(prob_vec_obs_all_total)=1:limit

library(abc)
param_df=data.frame(radius=radius_vec)
abc_results=abc(target=prob_vec_sharing2[1:limit], 
                param=param_df, 
                sumstat=prob_vec_obs_all, 
                tol=0.05, hcorr=T,method="neuralnet")

#Generate images
pdf("colon_patch_abc_results.pdf",useDingbats = F)
plot(abc_results, param=param_df)
dev.off()

pdf("colon_sharing_hist.pdf",height=5,width=8)
hist(round(colon_dists_vec),breaks=40,col='steelblue',xlab="Number of shared SNVs",main="",probability = T)
lines(d,lty='dashed')
abline(v=15,lwd=2,col='red')
dev.off()

pdf("crypt_distance_genetic_all_cutoff.pdf")
image(k, col=r,xlab="Physical distance (in crypts)",ylab="Number of shared SNVs",xlim=c(min(k$x),40),ylim=c(min(k$y),80))
abline(h=15,lwd=2,col='red')
dev.off()

