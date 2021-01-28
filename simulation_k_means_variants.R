
library(MASS)


# functions to perform k-means algorithm and RDC transformation
RDC_transform_s<- function (x, s="original", k=10, seed=19){
  k<- max(k, dim(x)[2])
  if(s=="original"){    s<-1/(ncol(x)+1)  }
  
  
  if(s==".3"){    s<-1/sqrt(ncol(x)+1)  }
  
  
  if(s==".5"){s<-3/(2*sqrt(ncol(x)+1))}
  
  set.seed(seed)
  
  x <- cbind(apply(as.matrix(x),2,function(u)ecdf(u)(u)),1)
  
  weights<- matrix(rnorm(ncol(x)*k),ncol(x))
  x <- s*x%*%weights
  x <- sin(x)
  sinx <- sin(x)
  s_weights<-s*weights
  result<- list(sinx=sinx, s_weights=s_weights)
  return(result)
}


RDC_transform_minsin<- function (x, s="original", k=10){
  
  if(s=="original"){    s<-1/ncol(x)  }
  if(s==".3"){    s<-1/sqrt(ncol(x))  }
  if(s==".5"){s<-3/(2*sqrt(ncol(x)))}
  
  
  k<- max(k, dim(x)[2])
  set.seed(2589)
  x <- cbind(apply(as.matrix(x),2,function(u)ecdf(u)(u)),1)
  x[,bin]<- xbin
  
  weights<- matrix(rnorm(ncol(x)*k),ncol(x))
  x <- s*x%*%weights
  s_weights<-s*weights
  result<- list(x=x, s_weights=s_weights)
  return(result)
}


kmeans_RDC<- function(x, s="original", k=10, seed=19, bin=0){
  
  x_RDC<- RDC_transform_s(x=x,s=s,k=k, seed=seed)$sinx
  kmeansres<- kmeans(x_RDC,2, nstart = 5)
  return(kmeansres)
}



kmeans_RDCminsin<- function(x, s="original", k=10, seed=19){
  
  x_RDC<- RDC_transform_minsin(x=x,s=s,k=k, seed=seed)$x
  kmeansres<- kmeans(x_RDC,2, nstart = 5)
  return(kmeansres)
}

kmeans_RDCmax<- function(x, s="original", k=10, seed=19){
  x_RDC<- RDC_transform_s(x=x,s=s,k=k, seed=seed)$sinx
  score<- c()
  for(i in 1:k){
    kmeansres<- kmeans(x_RDC[,i],2, nstart = 5)
    score[i]<- mean( kmeansres$withinss)/kmeansres$betweenss
  }
  best<- which(score==min(score))
  kmeansres<- kmeans(x_RDC[,best],2, nstart = 5)
  return(kmeansres)
  
  
}


kmeans_std<- function(x){
  
  std<-apply(x,2,function(u){(u-mean(u)/sd(u))})
  kmeansres<- kmeans(std,2, nstart = 5)
  return(kmeansres)
  
}

kmeans_plus<- function(x, plus=.5){
  n<- nrow(x)
  kmeansres<- kmeans(x,2, nstart = 5)
  score_raw<- kmeansres$betweenss/(kmeansres$betweenss+sum( kmeansres$withinss))
  
  if(score_raw<plus){
    score<- c()
    for(i in 1:ncol(x)){
      kmeansres<- kmeans(x[,i],2, nstart = 5)
      score[i]<- kmeansres$betweenss/(kmeansres$betweenss+sum( kmeansres$withinss))
    }
    
    weights<- matrix(score, n,ncol(x), byrow = T)
    x_PLUS<-   weights*x
    
    kmeansres_PLUS<- kmeans(x_PLUS,2, nstart = 5)
    scorePLUS<- mean(kmeansres$withinss)/kmeansres$betweenss
    
    if(scorePLUS<score_raw){return(kmeansres_PLUS)
    }else{return(kmeansres)}
  }else{return(kmeansres)}
  
}







kmeans_RDC_plus2<- function(x, s="original", k=10, plus=.5, seed=19){
  # This is another version of RDC plus where we use the explained variance, that did not work as well
  n<- dim(x)[1]
  # first run regular k means RDC
  # note that the plus component is not random
  x_RDC<- RDC_transform_s(x=x,s=s,k=k, seed=seed)$sinx
  kmeansres<- kmeans(x_RDC,2, nstart = 5)
  
  scoreRDC<- kmeansres$betweenss/(kmeansres$betweenss+sum( kmeansres$withinss))
  
  if(scoreRDC<plus){
    score<- c()
    for(i in 1:ncol(x)){
      kmeansres<- kmeans(x[,i],2, nstart = 5)
      score[i]<- kmeansres$betweenss/(kmeansres$betweenss+sum( kmeansres$withinss))
    }
    # assign the same weight to every score greater than or equal to 1
    # note that the within SS is exactly what the kmeans algorithm is estimating
    # we only use between to standardize it, but shouldnt change the results
    # as it is only multiplied by a constant (same when we take the mean of within SS)
    weights<- matrix(score, n,ncol(x), byrow = T)
    x_ecdf <- apply(as.matrix(x),2,function(u)ecdf(u)(u))
    x_PLUS<-   weights*x_ecdf
    
    kmeansres_PLUS<- kmeans(x_PLUS,2, nstart = 5)
    scorePLUS<- mean(kmeansres$withinss)/kmeansres$betweenss
    
    # we say if they are equal, return the multivariate split
    # so we use two different types of weights
    # and for now we assume that the sin step is unnecessary
    
    if(scorePLUS>scoreRDC){return(kmeansres_PLUS)
    }else{return(kmeansres)}
  }else{return(kmeansres)}
  
}


kmeans_RDC_plus<- function(x, s="original", k=10,  seed=19, plus=0){
  n<- dim(x)[1]
  # first run regular k means RDC
  # note that the plus component is not random
  x_RDC<- RDC_transform_s(x=x,s=s,k=k, seed=seed)$sinx
  kmeansres<- kmeans(x_RDC,2, nstart = 5)
  
  scoreRDC<- mean( kmeansres$withinss)/kmeansres$betweenss
  
  if(scoreRDC>plus){
    score<- c()
    for(i in 1:ncol(x)){
      kmeansres<- kmeans(x[,i],2, nstart = 5)
      score[i]<-mean( kmeansres$withinss)/kmeansres$betweenss
    }
    # assign the same weight to every score greater than or equal to 1
    # note that the within SS is exactly what the kmeans algorithm is estimating
    # we only use between to standardize it, but shouldnt change the results
    # as it is only multiplied by a constant (same when we take the mean of within SS)
    
    score[score>1]<-1
    weights<- matrix(1-score, n,ncol(x), byrow = T)
    x_ecdf <- apply(as.matrix(x),2,function(u)ecdf(u)(u))
    
    
    x_PLUS<-   weights*x_ecdf
    
    kmeansres_PLUS<- kmeans(x_PLUS,2, nstart = 5)
    scorePLUS<- mean(kmeansres$withinss)/kmeansres$betweenss
    
    # we say if they are equal, return the multivariate split
    # so we use two different types of weights
    # and for now we assume that the sin step is unnecessary
    
    if(scorePLUS<scoreRDC){return(kmeansres_PLUS)
    }else{return(kmeansres)}
  }else{return(kmeansres)}
  
}




kmeans_RDChalfmax<- function(x, s="original", k=10, seed=19){
  x_RDC<- RDC_transform_s(x=x,s=s,k=k, seed=seed)$sinx
  score<- c()
  for(i in 1:k){
    kmeansres<- kmeans(x_RDC[,i],2, nstart = 5)
    # we save the score, we want the between cluster variance to be large
    # and the within cluster variance to be small
    # so we want the score to be small
    score[i]<- mean(kmeansres$withinss)/kmeansres$betweenss
  }
  
  x_RDC[,order(score)]
  
  half_best<- round(k/2)
  
  kmeansres<- kmeans(x_RDC[,1:half_best],2, nstart = 5)
  return(kmeansres)
  
  
}

# functions to evaluate performance of 




performance<- function(kmeans, n=200){
  trueclusters2<-c(rep(2,n/2), rep(1,n/2))
  true_clusters<-c(rep(1,n/2), rep(2,n/2))
  max(sum(kmeans$cluster==true_clusters), sum(kmeans$cluster==trueclusters2))/n
}





# functions to generate data


normal_sample<-function(d=2,n=200,mua=1, mub=0, sda=1, sdb=1, rhoa=0.3, rhob=0.3, seed=346){
  # rho is the correlation 
  # a and b indicate the different groups
  # Parameters for bivariate normal distribution
  
  # change n such that it represents the group size
  n<-n/2
  mua<- rep(mua,d)[1:d]
  mub<- rep(mub,d)[1:d]
  
  Sigmaa<- matrix(rep(sda^2*rhoa, d^2),d,d)
  diag(Sigmaa)<- rep(sda^2,d)
  
  Sigmab<- matrix(rep(sdb^2*rhob, d^2),d,d)
  diag(Sigmab)<- rep(sdb^2,d)
  set.seed(seed)
  
  Xa <- mvrnorm(n=n, mu = mua, Sigma = Sigmaa ) 
  Xb <- mvrnorm(n=n, mu = mub, Sigma = Sigmab ) 
  
  x<- rbind(Xa,Xb)
  
  return(x)
}

normal_bin_noise_sample2<-function(d=2,n=200,mua=1, mub=0 , sda=1, sdb=1,p=.5, seed=346){

  Xa <- rnorm(n/2,mua,sda) 
  Xb <- rnorm(n/2,mub,sdb) 
  norm<- c(Xa,Xb)
  
  bin<- rbinom((d-1)*n, 1, p)
  
  x<- cbind(norm, matrix(bin, n, d-1))
  
  
  return(x)
}


normal_bin_noise_sample<-function(d=2,n=200,mua=1, mub=0 , sda=1, sdb=1,p=.5, seed=346){

  Xa <- rnorm(n/2,mua,sda) 
  Xb <- rnorm(n/2,mub,sdb) 
  norm<- c(Xa,Xb)
  
  x<- cbind(norm, c(rep(0,n/2),rep(1,n/2)))
  
  if(d>2){
    bin<- rbinom((d-2)*n, 1, p)
    x<- cbind(x, matrix(bin, n, d-2))
    
  }
  return(x)
}

normal_unif_noise_sample<-function(d=2,n=200,mua=1, mub=0, sda=1, sdb=1,unif1=2, unif2=5, seed=346){

  Xa <- rnorm(n/2,mua,sda) 
  Xb <- rnorm(n/2,mub,sdb) 
  norm<- c(Xa,Xb)
  
  
  unifs<- runif((d-1)*n, unif1, unif2)
  x<- cbind(norm, matrix(unifs, n, d-1))
  
  return(x)
}


# this one is with correlated noise
# more difficult to evaluate, as it is not straightforward which groups are best
bin_norm_noise_cor_sample<-function(d=2,n=200,mu=1, sd=1, rho=0.3, p=.5, seed=346){
  

  
  bin<- c(rep(0,n*p), rep(1,n*p))
  # change d such that it indicates the total number of variables including the normal 
  d<- d-1
  set.seed(seed)
  
  if(d==1){
    norm<- rnorm(n, mu, sd)
  }else{
    mu<- rep(mu,d)
    
    Sigma<- matrix(rep(rho,d^2),d,d)
    diag(Sigma)<- sd
    
    norm <- mvrnorm(n=n, mu = mu, Sigma = Sigma ) 
  }
  
  x<- cbind(bin,norm)
  
  return(x)
}



bin_norm_noise_sample<-function(d=2,n=200,mu=5, sd=2, rho=0.3, p=.5, seed=346){
  

  
  bin<- c(rep(0,n*p), rep(1,n*p))
  # change d such that it indicates the total number of variables including the normal 
  d<- d-1
  set.seed(seed)
  
  
  norm<- matrix(rnorm(n*d, mu, sd), n, d)
  
  
  x<- cbind(bin,norm)
  
  return(x)
}



bin_exp_noise_sample<-function(d=2,n=200,expa=1/14,p=.5, seed=346){

  
  bin<- c(rep(0,n*p), rep(1,n*p))
  # change d such that it indicates the total number of variables including the normal 
  d<- d-1
  set.seed(seed)
  
  
  exp<- matrix(rexp(n*d, expa), n, d)
  
  
  x<- cbind(bin,exp)
  
  return(x)
}







poiss_sample<-function(d=2,n=200,poissa=2, poissb=4 ,seed=346){

  
  set.seed(seed)
  
  Xa <- matrix(rpois(n=n/2*d, poissa), n/2,d) # from MASS package
  Xb <- matrix(rpois(n=n/2*d, poissb), n/2,d) # from MASS package
  
  x<- rbind(Xa,Xb)
  
  return(x)
}






exp_sample<-function(d=2,n=200,expa=1/14, expb=1/7 , delta=0, seed=346){

  
  set.seed(seed)
  
  Xa <- matrix(rexp(n=n/2*d, expa), n/2,d) # from MASS package
  Xb <- matrix(rexp(n=n/2*d, expb)+delta, n/2,d) # from MASS package
  
  x<- rbind(Xa,Xb)
  
  return(x)
}



exp_bin_noise_sample<-function(d=2,n=200,expa=1/14, expb=1/7 , p=.5, delta=0, seed=346){

  set.seed(seed)
  
  Xa <- rexp(n=n/2, expa) # from MASS package
  Xb <- rexp(n=n/2, expb)+delta # from MASS package
  
  exp<- c(Xa,Xb)
  x<- cbind(exp, c(rep(0,n/2),rep(1,n/2)))
  
  
  if(d>2){
    bin<- rbinom((d-2)*n, 1, p)
    x<- cbind(x, matrix(bin, n, d-2))
    
  }
  return(x)
}




exp_unif_noise_sample<-function(d=2,n=200,expa=1/14, expb=1/7 , unif1=2, unif2=5, delta=0, seed=346){
  # scheve verdeling  (niet normaal)
  # groepen zitten in de continue variabele
  # noise is ook continu (uniform)
  
  set.seed(seed)
  
  Xa <- rexp(n=n/2, expa) # from MASS package
  Xb <- rexp(n=n/2, expb)+delta # from MASS package
  
  exp<- c(Xa,Xb)
  
  unifs<- runif((d-1)*n, unif1, unif2)
  x<- cbind(exp, matrix(unifs, n, d-1))
  return(x)
}



K_means_results<- function(x, seed=19, plus=0){
  n<- nrow(x)
  kmeans_results<- matrix(NA, 1,10)
  
  kmeans_results[1]<- performance(kmeans_RDC(x, seed=seed), n=n)
  
  kmeans_results[2]<- performance(kmeans_RDC(x, s=".3", seed=seed), n=n)
  
  kmeans_results[3]<- performance(kmeans_RDC(x,s=".5", seed=seed), n=n)
  
  kmeans_results[4]<-performance(kmeans_RDCmax(x, seed=seed), n=n)
  
  kmeans_results[5]<-performance(kmeans_RDCmax(x, s=".5", seed=seed), n=n)
  
  kmeans_results[6]<-performance(kmeans_RDC_plus2(x, seed=seed), n=n)
  
  kmeans_results[7]<- performance(kmeans_std(x), n=n)
  
  kmeans_results[8]<-performance(kmeans(x,2, nstart=5), n=n)
  
  kmeans_results[9]<-performance(kmeans_plus(x,plus=plus), n=n)
  
  kmeans_results[10]<-performance(kmeans_RDC_plus(x), n=n)
  
  colnames(kmeans_results)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus2", "standardized k-means", "raw k-means", "kmeans plus", "RDC+")
  
  kmeans_results
  
}


K_means_simulation_norm<- function(n=200,mua=1, mub=0, sda=1, sdb=1, rhoa=0.3, rhob=0.3, seed=346, sims=100, scores=F){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- normal_sample(d=d, mua=mua, mub=mub, rhoa=rhoa, rhob=rhob, n=n, sda=sda, sdb=sdb, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}



K_means_simulation_norm_unif_noise<- function(n=200, mua=1, mub=0, sda=1, sdb=1,unif1=2, unif2=5, seed=346, sims=100, k=10, check_linear=F, delta=0, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- normal_unif_noise_sample(d=d, mua=mua, mub=mub, sda=sda, n=n, sdb=sdb, unif1=unif1, unif2=unif2, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}





K_means_simulation_norm_bin_noise<- function(n=200,mua=1,  mub=0, sda=1, sdb=1, p=0.5, seed=346, sims=100, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- normal_bin_noise_sample(d=d, mua=mua, mub=mub, sda=sda, n=n, sdb=sdb, p=p, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}

################################################################################


K_means_simulation_bin_norm_noise<- function(n=200,mu=1,  sd=1, rho=0.3, p=0.5, seed=346, sims=100, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- bin_norm_noise_sample(d=d, mu=mu,rho=rho, n=n, sd=sd, p=p, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}


K_means_simulation_bin_exp_noise<- function(n=200, expa=1/14, p=.5, seed=346, sims=100, k=10, check_linear=F, delta=0, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- bin_exp_noise_sample(d=d, expa=expa, p=p, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}

################################################################################



K_means_simulation_exp<- function(n=200,expa=1/14, expb=1/7 , delta=0,  seed=346, sims=100, k=10, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- exp_sample(d=d, expa=expa, expb=expb, n=n, delta=delta, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}


K_means_simulation_exp_unif_noise<-function(n=200,expa=1/14, expb=1/7 , unif1=2, unif2=5,
                                            seed=346, sims=100, k=10, check_linear=F, delta=0, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- exp_unif_noise_sample(d=d, expa=expa,  expb=expb, unif1=unif1, unif2=unif2,delta=delta, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}







K_means_simulation_exp_bin_noise<- function(n=200, expa=1/14, expb=1/7 , p=.3 ,
                                            seed=346, sims=100, k=10, check_linear=F, delta=0, scores=FALSE){
  
  k_means_results<- array(0,c(19,sims,10))
  
  for(d in 2:20){  
    for(j in 1: sims){
      x<- exp_bin_noise_sample(d=d, expa=expa,  expb=expb, p=p, delta=delta, seed=seed+j+d)
      if(scores==FALSE){
        k_means_results[d-1,j,]<- K_means_results(x, seed=seed+j+d)
      }else{k_means_results[d-1,j,]<- K_means_scores(x)}
      
    }  
  }
  
  result<- apply(k_means_results,c(1,3),median)
  colnames(result)<- c("RDC", "RDC s=0.3", "RDC s=0.5",  "RDC max", "RDC max s=0.5", "RDC plus", "standardized k-means", "raw k-means", "kmeans plus", "RDC+ old")
  result
}


###################################################################################





###########################################################################################################################

# Function to repeat above for multiple thresholds for RDC

test_plus_threshold_old<- function(n=200, sims=250, plus= seq(.2,2,.05), seed=346){
  
  # use the defaults for all the data generation 
  
  threshold_results_RDC<- array(NA, c(19,8, length(plus)))
  
  
  for(i in 1:length(plus)){
    for(d in 2:20){
      
      x<- normal_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,1,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- normal_bin_noise_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,2,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- normal_unif_noise_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,3,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- exp_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,4,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- exp_bin_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,5,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- exp_uniform_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,6,i]<-performance(kmeans_RDC_plus(x,plus=plus[i], seed=seed+d+i), n=n)
      
      x<- bin_norm_noise_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,7,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      x<- bin_exp_noise_sample(n=n, d=d, seed=seed+d+i)
      threshold_results_RDC[d-1,8,i]<-performance(kmeans_RDC_plus(x, plus=plus[i], seed=seed+d+i), n=n)
      
      
    }
  }
  threshold_results_RDC
}



















###########################################################################################################################

# Let's first see which threshold is fit

# Try to find good threshold
threshold200<- test_plus_threshold_old(sims=1000, plus=seq(0,2,.05))
threshold1000<- test_plus_threshold_old(n=1000,sims=1000, plus=seq(0,2,.05))

plus<- seq(0,2,.05)
scenarios<- c("norm_overlap","norm_bin_nois","norm_unif_noise","exp_overlap","exp_bin_noise","exp_unif_noise","bin_norm","bin_exp")
rownames(threshold1000)<- rownames(threshold200)<- 2:20



# create a heatmap of the results, with the performance as the colors
# varying over the threshold and the number of data points

# we are primarily interested in the mixed scenarios: these are scenario 2,5,7 and 8
# we can mark these scenarios to accentuate them: 
scenarios<- c("norm_overlap","!norm_bin_nois","norm_unif_noise","exp_overlap","!exp_bin_noise","exp_unif_noise","!bin_norm","!bin_exp")

for(i in 1:8){
  dat<-threshold200[,i,]
  colnames(dat)<- plus
  rownames(dat)<- as.character(2:20)
  heatmap(dat, Rowv=NA, Colv=NA, col=heat.colors(11), main=paste(c("n=200", scenarios[i]), collapse = " "))
  legend(xy.coords(.8,1.04) , legend=seq(0.5,1,.05),fill=heat.colors(11))
}
for(i in 1:8){
  dat<-threshold1000[,i,]
  colnames(dat)<- plus
  rownames(dat)<- as.character(2:20)
  heatmap(dat, Rowv=NA, Colv=NA, col=heat.colors(11), main=paste(c("n=1000", scenarios[i]), collapse = " "))
  legend(xy.coords(.8,1.04) , legend=seq(0.5,1,.05),fill=heat.colors(11))
}


###########################################################################################################################

###########################################################################################################################



# Results: compare methods 

norm_overlap<- K_means_simulation_norm(sims=500)


norm_bin_noise<- K_means_simulation_norm_bin_noise(sims=500)


norm_unif_noise<- K_means_simulation_norm_unif_noise(sims=500)


exp_bin_noise<- K_means_simulation_exp_bin_noise(sims=500)


exp_unif_noise<- K_means_simulation_exp_unif_noise(sims=500)


exp_overlap<- K_means_simulation_exp(sims=500)


bin_norm<- K_means_simulation_bin_norm_noise(sims=500) 


bin_exp<- K_means_simulation_bin_exp_noise(sims=500) 

###########################################################################################################################

matplot(norm_bin_noise[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1 )


matplot(exp_bin_noise[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)



matplot(bin_norm[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)


matplot(bin_exp[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)




# only for continuous scenarios

matplot(norm_overlap[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend(xy.coords(.8,1.04), legend = trans_names,col=1:7, pch=1)


matplot(norm_unif_noise[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)


matplot(exp_overlap[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)


matplot(exp_unif_noise[,trans], type = "b",pch=1,col = 1:7, xlab = "Number of variables", ylab = "Proportion of correct classifications")
legend("topright", legend = trans_names,col=1:7, pch=1)
