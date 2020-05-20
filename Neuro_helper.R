Simu_data <- function(x,y,
                      bandwidth,
                      pic_x,pic_y,
                      width,
                      amplitude,
                      noisesd){
  xx <- pic_x - x
  yy <- pic_y - y
  
  absorb <- (pnorm(xx+width,sd = bandwidth) - pnorm(xx-width,sd = bandwidth)) * 
    (pnorm(yy+width,sd = bandwidth) - pnorm(yy-width,sd = bandwidth))
  
  rnorm(1, mean = amplitude * absorb, sd = noisesd)
}

pred_data <- function(x,y,
                      bandwidth,
                      pic_x,pic_y,
                      width,
                      amplitude
                      ){
  xx <- pic_x - x
  yy <- pic_y - y
  
  absorb <- (pnorm(xx+width,sd = bandwidth) - pnorm(xx-width,sd = bandwidth)) * 
    (pnorm(yy+width,sd = bandwidth) - pnorm(yy-width,sd = bandwidth))
  
  amplitude * absorb
}

log_Lik_single <- function(x,y,
                           bandwidth,
                           pic_x,pic_y,
                           width,
                           amplitude,
                           noisesd,obs){
  
  xx <- pic_x - x
  yy <- pic_y - y
  
  absorb <- (pnorm(xx+width,sd = bandwidth) - pnorm(xx-width,sd = bandwidth)) * 
    (pnorm(yy+width,sd = bandwidth) - pnorm(yy-width,sd = bandwidth))
  
  dnorm(obs, mean = amplitude * absorb, sd = noisesd,log = T)
  
}


Simu_data_batch <- function(x,y,
                            bandwidth,
                            pic_pos,
                            width,
                            amplitude,
                            noisesd){
  
  n <- nrow(pic_pos)
  sapply(1:n,function(i,x,y,
                      bandwidth,
                      pic_pos,
                      width,
                      amplitude,
                      noisesd){
    Simu_data(x,y,
              bandwidth,
              pic_pos[i,1],pic_pos[i,2],
              width,
              amplitude,
              noisesd)
  },x,y,
  bandwidth,
  pic_pos,
  width,
  amplitude,
  noisesd)
}

pred_batch <- function(x,y,
                            bandwidth,
                            pic_pos,
                            width,
                            amplitude
                            ){
  
  n <- nrow(pic_pos)
  sapply(1:n,function(i,x,y,
                      bandwidth,
                      pic_pos,
                      width,
                      amplitude,
                      noisesd){
    pred_data(x,y,
              bandwidth,
              pic_pos[i,1],pic_pos[i,2],
              width,
              amplitude
              )
  },x,y,
  bandwidth,
  pic_pos,
  width,
  amplitude
  )
}

log_Lik_batch <- function(obs_mat,x,y,
                            bandwidth,
                            width,
                            amplitude,
                            noisesd){
  
  n <- nrow(obs_mat)
  sum(sapply(1:n,function(i,x,y,
                      bandwidth,
                      obs_mat,
                      width,
                      amplitude,
                      noisesd){
    log_Lik_single(x,y,
              bandwidth,
              obs_mat[i,1],obs_mat[i,2],
              width,
              amplitude,
              noisesd,obs_mat[i,3])
  },x,y,
  bandwidth,
  obs_mat,
  width,
  amplitude,
  noisesd))
}

prior_xy_unif <- function(x,y,log = TRUE){
  ifelse(log,dunif(x,log = T) + dunif(y,log = T),
         dunif(x,log = F) * dunif(y,log = F))
  
  
}

Main_sampler <- function(obs_mat, width, alpha, beta, prior_xy = prior_xy_unif,
                         n_iter=1000,burn_in=100,
                         thin_by = 1,prop_sd_posi = 0.01,
                         prop_sd_amplitude = .5,
                         prop_sd_bandwidth = 0.015){
  
  require(coda)
  n_save <- floor( n_iter/thin_by)
  noisesd_mcmc <- mcmc(matrix(NA,nrow = n_save))
  amplitude_mcmc <- mcmc(matrix(NA,nrow = n_save))
  posi_mcmc <- mcmc(matrix(NA,nrow = n_save,ncol = 2))
  bandwidth_mcmc <- mcmc(matrix(NA,nrow = n_save))
  log_lik_mcmc <- mcmc(matrix(NA,nrow = n_save))
  
  noisesd_curr <- sqrt( 1/rgamma(1,alpha,scale = beta)) # noise sd
  amplitude_curr <- rexp(1,0.1) # initial value for amplitude
  posi_curr <- runif(2) # initial position 
  bandwidth_curr <- rexp(1,1) # initial bandwidth
  
  log_posterior_curr <- log_Lik_batch(obs_mat,posi_curr[1],posi_curr[2],
                                      bandwidth_curr,
                                      width,
                                      amplitude_curr,
                                      noisesd_curr) + 
                        prior_xy(x=posi_curr[1],y=posi_curr[2]) + 
                        dgamma(1/noisesd_curr^2,alpha,beta,log = T)
    
  k <- 1 # how many sample have been saved
  for(i in 1:(n_iter+burn_in)) {
    ## update noisesd using conjugate
    mean_pred_curr <- pred_batch(posi_curr[1],posi_curr[2],
                           bandwidth_curr,
                           obs_mat[,1:2],
                           width,
                           amplitude_curr)
    SSE <- sum((mean_pred_curr-obs_mat[,3])^2)
    
    noisesd_curr <- sqrt( 1/(rgamma(1,alpha+length(mean_pred_curr)/2,
                              beta + SSE/2)))
    log_posterior_curr <- log_Lik_batch(obs_mat,posi_curr[1],posi_curr[2],
                                        bandwidth_curr,
                                        width,
                                        amplitude_curr,
                                        noisesd_curr) + 
      prior_xy(x=posi_curr[1],y=posi_curr[2])+ 
      dgamma(1/noisesd_curr^2,alpha,beta,log = T)
    
    ## update amplitude using MH
    amplitude_prop <- rnorm(1,amplitude_curr,prop_sd_amplitude)
    if(amplitude_prop>0){
      log_posterior_prop <-  log_Lik_batch(obs_mat,posi_curr[1],posi_curr[2],
                                 bandwidth_curr,
                                 width,
                                 amplitude_prop,
                                 noisesd_curr) + 
                            prior_xy(x=posi_curr[1],y=posi_curr[2])+ 
                            dgamma(1/noisesd_curr^2,alpha, beta,log = T)
      
      if(log(runif(1))<log_posterior_prop-log_posterior_curr){
        amplitude_curr <- amplitude_prop
        log_posterior_curr <- log_posterior_prop
      }
    }
    
    ## updating bandwidth using MH
    bandwidth_prop <- rnorm(1,bandwidth_curr,prop_sd_bandwidth)
    if(bandwidth_prop>0){
      log_posterior_prop <-  log_Lik_batch(obs_mat,posi_curr[1],posi_curr[2],
                                           bandwidth_prop,
                                           width,
                                           amplitude_curr,
                                           noisesd_curr) + 
        prior_xy(x=posi_curr[1],y=posi_curr[2])+ 
        dgamma(1/noisesd_curr^2,alpha,beta,log = T)
      
      if(log(runif(1))<log_posterior_prop-log_posterior_curr){
        bandwidth_curr <- bandwidth_prop
        log_posterior_curr <- log_posterior_prop
      }
    }
    
    ## updating position
    posi_prop <- rnorm(2,posi_curr,prop_sd_posi)
    if(all(posi_prop>=0 & posi_prop<=1)){
      log_posterior_prop <-  log_Lik_batch(obs_mat,posi_prop[1],posi_prop[2],
                                           bandwidth_curr,
                                           width,
                                           amplitude_curr,
                                           noisesd_curr) + 
        prior_xy(x=posi_prop[1],y=posi_prop[2])+ 
        dgamma(1/noisesd_curr^2,alpha,beta,log = T)
      
      if(log(runif(1))<log_posterior_prop-log_posterior_curr){
        posi_curr <- posi_prop
        log_posterior_curr <- log_posterior_prop
      }
    }
    
    if( i>burn_in &  (i-burn_in) %% thin_by==0){
      cat(k/n_save,"\n")
      noisesd_mcmc[k,] <- noisesd_curr
      amplitude_mcmc[k,] <- amplitude_curr
      posi_mcmc[k,] <- posi_curr
      bandwidth_mcmc[k,] <- bandwidth_curr
      log_lik_mcmc[k,] <- log_Lik_batch(obs_mat,posi_curr[1],posi_curr[2],
                                        bandwidth_curr,
                                        width,
                                        amplitude_curr,
                                        noisesd_curr)
      k <- k+1
      
    }
  }
  
  return(list(noise = noisesd_mcmc,
              amplitude = amplitude_mcmc,
              bandwidth = bandwidth_mcmc,
              posi = posi_mcmc,
              log_lik = log_lik_mcmc))
  
}






