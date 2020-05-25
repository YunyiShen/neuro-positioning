source("./R/Neuro_helper.R")

bandwidth <- 0.1
width <- .2

posi <- c(.3,.5)

set.seed(786)



experimentxy <- as.matrix( expand.grid(seq(0,1,width),seq(0,1,width)))

amplitude <- 10
noisesd <- 0.3

n_rep <- 5

simu_list <- lapply(1:n_rep,function(i,posi,
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

Mean_response <- lapply(simu_list,function(w){
  matrix(w,nrow = sqrt(length(w)))
})

Mean_response <- Reduce("+",Mean_response)
Mean_response <- Mean_response/n_rep
Mean_response <- Mean_response/sum(Mean_response)
image(Mean_response)

obs_list <- lapply(simu_list,cbind,experimentxy)
obs_mat <- Reduce(rbind,obs_list)
obs_mat <- obs_mat[,c(2,3,1)]

first_test <- Main_sampler(obs_mat, width, alpha = 1, beta = 1, prior_xy = prior_xy_unif,
                           n_iter=5000,burn_in=3000,
                           thin_by = 10,prop_sd_posi = 0.01, # change this to 0.05 if it does not work
                           prop_sd_amplitude = .5,
                           prop_sd_bandwidth = 0.01)

dev.off()
par(mfrow = c(1,2))

image(Mean_response,main = "mean response")
points(as.matrix(first_test$posi))
points(x = posi[1],y = posi[2],col = "red",pch = 18)

densplot(first_test$bandwidth,main = "bandwidth",show.obs = T,xlab = "")
abline(v = bandwidth,col = "red")

