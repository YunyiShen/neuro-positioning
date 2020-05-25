if(!require(svMisc)) install.packages("svMisc")
if(!require(coda)) install.packages("coda")
source("./R/Neuro_helper.R")



bandwidth_list <- c(0.1,0.15,0.2)
width_list <- c(0.1,0.125,0.2,0.25,0.5)
#bandwidth_list <- c(0.2)
#width_list <- c(0.2)

posi <- c(.3,.5)
amplitude <- 10
noisesd <- 0.3

set.seed(786)

n_rep <- 5
n_data_set <- 150
#n_data_set <- 2

cat("\n\nstart the simulation\n number of perception bandwidth:",
    length(bandwidth_list),
    "\n number of picture width:",
    length(width_list),
    "\n number of data sets for each setting:",
    n_data_set,"\n\n")

for(bandwidth in bandwidth_list){
  for(width in width_list){
    a_post <- gamma_post <- y_c_post <- x_c_post <- sigma_post <- 
      matrix(NA,nrow = n_data_set,ncol = 500)
    experimentxy <- as.matrix( expand.grid(seq(0,1,width),seq(0,1,width)))
    for(i in 1:n_data_set){
      cat("perception bandwidth =",bandwidth," picture width =",width," dataset number:",i,"\n\n")
      simu_list <- lapply(1:n_rep,function(foo,posi,
                                           bandwidth,
                                           experimentxy,width,
                                           amplitude,noisesd){
        Simu_data_batch(posi[1],posi[2],
                        bandwidth,
                        experimentxy,width,
                        amplitude,noisesd)
      },posi,
      bandwidth,
      experimentxy,width,
      amplitude,noisesd
      )
      
      obs_list <- lapply(simu_list,cbind,experimentxy)
      obs_mat <- Reduce(rbind,obs_list)
      obs_mat <- obs_mat[,c(2,3,1)]
      
      first_test <- Main_sampler(obs_mat, width, alpha = 1, beta = 1, prior_xy = prior_xy_unif,
                                 n_iter=5000,burn_in=3000,
                                 thin_by = 10,prop_sd_posi = 0.01, # change this to 0.05 if it does not work
                                 prop_sd_amplitude = .5,
                                 prop_sd_bandwidth = 0.01)
      
      a_post[i,] <- first_test$amplitude
      gamma_post[i,] <- first_test$bandwidth
      x_c_post[i,] <- first_test$posi[,1]
      y_c_post[i,] <- first_test$posi[,2]
      sigma_post[i,] <- first_test$noise
      
      write.csv(a_post,paste0("./Simulations/Big_simu/bw_",bandwidth,"/w_",width,"/a.csv"))
      write.csv(gamma_post,paste0("./Simulations/Big_simu/bw_",bandwidth,"/w_",width,"/gamma.csv"))
      write.csv(x_c_post,paste0("./Simulations/Big_simu/bw_",bandwidth,"/w_",width,"/x_c.csv"))
      write.csv(y_c_post,paste0("./Simulations/Big_simu/bw_",bandwidth,"/w_",width,"/y_c.csv"))
      write.csv(sigma_post,paste0("./Simulations/Big_simu/bw_",bandwidth,"/w_",width,"/sigma.csv"))
      
    }
    
  }
  
}
