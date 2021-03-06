source("./R/Neuro_helper.R")
posi <- c(0.3,0.5)
bandwidth <- 0.2
width <- 0.1

experimentxy <- as.matrix( expand.grid(seq(0,1,width),seq(0,1,width)))

amplitude <- 8
noisesd <- 0.3

large_simu <- Simu_data_batch(posi[1],posi[2],
                              bandwidth,
                              experimentxy,width,
                              amplitude,noisesd)

image(matrix(large_simu,sqrt(length(large_simu))))

obs_mat <- cbind(experimentxy,large_simu)


log_Lik_batch(obs_mat,posi[1],posi[2],
              bandwidth,
              width,
              amplitude,
              noisesd)

pred <- pred_batch(posi[1],posi[2],
           bandwidth,
           experimentxy,
           width,
           amplitude
           )
first_test <- Main_sampler(obs_mat, width, alpha = 1, beta = 1, prior_xy = prior_xy_unif,
                         n_iter=5000,burn_in=3000,
                         thin_by = 10,prop_sd_posi = 0.01, # change this to 0.05 if it does not work
                         prop_sd_amplitude = .5,
                         prop_sd_bandwidth = 0.01)

image(matrix(large_simu,sqrt(length(large_simu))))
points(as.matrix(first_test$posi))
points(x = 0.3,y = 0.5,col = "blue")

plot(first_test$amplitude)
plot(first_test$bandwidth)
plot(first_test$log_lik)
