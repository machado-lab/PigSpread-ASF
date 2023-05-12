## PEDV CALIBRATION ####

###############
#### Setup ####
###############

##Libraries
library(mc2d)
library(lubridate)
library(dplyr)
library(SimInf)
library(doSNOW)
library(doParallel)
library(foreach)

##Load data
load("~/PED_epidemic_data.RData")
load("~/population_and_movement_data.RData")

#################################
#### Define model parameters ####
#################################

biosecurity <- id_nodes$farm_type
bio_types <- unique(biosecurity)

tE_time = 5 #median days farms may be exposed

farm_detection_ratio = 10 #half of the average detection time which is 20 days
efetive_surveillance <- rep(0, nrow(matrix_1))

#Time step and start date
tSim = 140 # number of simulated time steps
time_step = time_step # unit of time from the time steps
date_start = as.Date("2020-01-13") # start date model
time_period = seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by="days")

##Farm nodes and movement/truck networks
nodes = id_nodes
movement_dataframe = db_net # database with network
movement_matrix = network_pig_movement # list of movement matrix by each time step
truck_feed_matrix = network_truck_feed
truck_pig_matrix = network_truck_pig
truck_market_matrix = network_truck_market
truck_loading_matrix = network_truck_loading

error <- c(20, # mu_total_cases_sow
           1, # mu_median_cases_sow
           10, # mu_max_cases_sow
           
           25, # mu_total_cases_fin
           1, # mu_median_cases_fin
           10, # mu_max_cases_fin
           
           20, # mu_total_cases_nur
           1, # mu_median_cases_nur
           10, # mu_max_cases_nur
           
           10000 #mu_mean_dist (in meters)
)

##Spatial component
#cut-off values are in meters!
sp_cutoff <- c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000)

dmx_list <- list()

for (i in 1:length(sp_cutoff)) {
  
  dmx <- distance_matrix(x = id_nodes$long,
                         y = id_nodes$lat,
                         cutoff = sp_cutoff[i],
                         min_dist = 100)
  
  dmx <- mxt_pop/(dmx^2)
  dmx[dmx == Inf] <- 0
  diag(dmx) <- 0
  
  dmx_list[(length(dmx_list) + 1)] <- dmx
}

number_accepted <- NULL
beta_accepted <- NULL

###############
#### Model ####
###############

limit <- data.frame(limit_1 = seq(1, 180001, by = 20000),
                    limit_2 = seq(20000, 200000, by = 20000))

for(ru in 1:10){

cores <- 50
cl <- makeSOCKcluster(cores, outfile = "~/mylog.txt")
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=20000, style=2)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

ss <- foreach(w = limit$limit_1[ru]:limit$limit_2[ru],
              .options.snow=opts,
              .packages=c( "tidyverse", "Matrix",
                           "doParallel", "foreach", "sp", "mc2d", "lubridate", "dplyr", "SimInf"), 
              .verbose = T,
              .inorder = TRUE,
              .errorhandling = "pass"
              
) %dopar% {
  
  cat(w)

    date_start <- sample(time_period, 1)
    
    if(date_start <= as.Date("2020-08-14")){
      
      time_start <- date_start + seq(0, (tSim-1), by = 1) #creating date for each time step
    }else{
      time_1 <- as.numeric(max(time_period) - date_start)
      time_2 <- abs(tSim - time_1)
      time_start <- c(date_start + seq(0, time_1, by = 1), min(time_period) + seq(0, time_2 - 2, by = 1))
    }
    
    time_start <- as.numeric(yday(time_start)) #formatting
    
    E <- rep(x = 0, times = nrow(matrix_1))
    I <- rep(x = 0, times = nrow(matrix_1))
    D <- rep(x = 0, times = nrow(matrix_1))
    
    ##Create matrix to store farm status for each time step 
    M_Sim_E = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1)) #store exposed 
    M_Sim_I = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1)) #store infected
    M_Sim_D = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1)) #store removed
    
    ##Seed infection at random farm
    seed_infect <- sample(1:2294, size = 1)
    I[seed_infect] <- 1 #assign this farm as infected
    
    dmx_num <- sample(1:13, 1)
    gravity_matrix = dmx_list[[dmx_num]]
    
    ##Just storing information about which farms were infected and the time they will be infected for
    tE <- E 
    tI <- I
    tD <- D
    
    ##We add status information to the matrix at col 1 (which is t0)
    M_Sim_E[, 1] = tE
    M_Sim_I[, 1] = tI
    M_Sim_D[, 1] = tD 
    
    ##Create data frame to store the amount of transmission through each route
    trans_routes <- data.frame()
    
    ##Indicates which farm is infected - this will then be used to add on top of to create a vector of the number of days a farm has been infected
    time_infected <- as.numeric(I>0)
    time_exposed <- as.numeric(E>0)
    time_detected <- as.numeric(D>0)
    
    bio_sow <- 0
    bio_nur <- 0
    bio_fin <- 0
    bio_gilt <- 0
    bio_boar <- 0
    
    if(length(number_accepted) < 10){
      
      surveil_sow <- runif(1, 0.85, 0.95)
      surveil_nur <- runif(1, 0.15, 0.2)
      surveil_fin <- runif(1, 0.5, 0.6)
      surveil_gilt <- runif(1, 0, 0.5)
      surveil_boar <- runif(1, 0, 0.5)
      
    }else{
      
      surveil_sow <- rnorm(1, mean(surveil_values_sow), sd(surveil_values_sow))
      surveil_nur <- rnorm(1, mean(surveil_values_nur), sd(surveil_values_nur))
      surveil_fin <- rnorm(1, mean(surveil_values_fin), sd(surveil_values_fin))
      surveil_gilt <- rnorm(1, mean(surveil_values_gilt), sd(surveil_values_gilt))
      surveil_boar <- rnorm(1, mean(surveil_values_boar), sd(surveil_values_boar))
      
    }
    
    bio_values <- c(bio_fin, bio_nur, bio_sow, bio_gilt, bio_boar)
    
    for(i in 1:length(bio_types)){
      biosecurity[which(biosecurity == bio_types[i])] <- bio_values[i]
    }
    
    biosecurity <- as.numeric(biosecurity)
    
    efetive_surveillance[farm_sow] <- surveil_sow
    efetive_surveillance[farm_nursery] <- surveil_nur
    efetive_surveillance[farm_finisher] <- surveil_fin
    efetive_surveillance[farm_gilt] <- surveil_gilt
    efetive_surveillance[farm_others] <- surveil_boar
    
    beta_values1_E <- NULL
    beta_values1_I <- NULL
    net2 <- NULL
    net_D2 <- NULL
    
    epi_t = 1
    
    ##Start of for loop for time   
    for (i in time_start[2:length(time_start)]) { #
      
      epi_t <- epi_t + 1
      
      #Pig movements in time step i
      movement_matrix_timestep <- movement_matrix[[i]]
      
      #Truck movements in time step i
      truck_feed_matrix_timestep <- truck_feed_matrix[[i]] 
      truck_pig_matrix_timestep <- truck_pig_matrix[[i]] 
      truck_market_matrix_timestep <- truck_market_matrix[[i]] 
      truck_loading_matrix_timestep <- truck_loading_matrix[[i]] 
      
      #Identifying which farms are in which compartment - binary results so 1 = they are in that compartment. 
      tE = M_Sim_E[, epi_t - 1]
      tI = M_Sim_I[, epi_t - 1]
      tD = M_Sim_D[, epi_t - 1]
      
      E = as.numeric(tE > 0) 
      I = as.numeric(tI > 0)
      D = as.numeric(tD > 0)
      S = as.numeric(!E & !I & !D) 
      
      #Rate of detection of infected farms
      ratedetection <- (I) * ((0.95/(1 + exp(-0.5*(time_infected - farm_detection_ratio)))) * efetive_surveillance)
      
      #This is the time they have been infected/exposed/removed for 
      time_infected <- time_infected + as.numeric(tI>0)
      time_infected <- time_infected * I
      
      ratedetection <- ratedetection * (M_Sim_I[, epi_t - 1] >= 1)
      # here decided threshold to detect the farm
      threshold_detection <- runif(length(ratedetection), 0.5, 1)
      
      # identified new farm with outbreak
      detected_farm <- as.numeric(ratedetection>threshold_detection)
      
      time_exposed <- time_exposed + as.numeric(tE > 0)
      time_exposed <- time_exposed * E
      
      time_detected <- time_detected + as.numeric(tD > 0)
      time_detected <- time_detected * D
      
      #This indicates farms which the infected farm may have move to in this time step
      Pvector <- (I %*% movement_matrix_timestep)[1,]
      Pvector_E <- (E %*% movement_matrix_timestep)[1,]
      Pvector_D <- (D %*% movement_matrix_timestep)[1,]
      truck_feed_vector <- (I %*% truck_feed_matrix_timestep)[1,]
      truck_feed_vector_D <- (D %*% truck_feed_matrix_timestep)[1,]
      truck_pig_vector <- (I %*% truck_pig_matrix_timestep)[1,] 
      truck_pig_vector_D <- (D %*% truck_pig_matrix_timestep)[1,]
      truck_market_vector <- (I %*% truck_market_matrix_timestep)[1,]
      truck_market_vector_D <- (D %*% truck_market_matrix_timestep)[1,]
      truck_loading_vector <- (I %*% truck_loading_matrix_timestep)[1,]
      truck_loading_vector_D <- (D %*% truck_loading_matrix_timestep)[1,]
      
      spatial_vector <- I %*% gravity_matrix
      spatial_vector_D <- D %*% gravity_matrix
      
      #Draw betas from a distribution
      
      if(length(number_accepted) < 10){
        
        beta_net_dist <- rpert((length(Pvector) + length(Pvector_D)), min = 1.4, mode = 1.6, max = 1.8, shape = 4)
        
        beta_net_dist_E <- rpert(length(Pvector_E), min = 1.4, mode = 1.6, max = 1.8, shape = 4)
        
        beta_truck_feed_dist <- rpert((length(truck_feed_vector) + length(truck_feed_vector_D)), min = 0.000024, mode = 0.00003, max = 0.000034, shape = 4)
        
        beta_truck_pig_dist <- rpert((length(truck_pig_vector) + length(truck_pig_vector_D)), min = 0.00035, mode = 0.0006, max = 0.00085, shape = 4)
        
        beta_truck_market_dist <- rpert((length(truck_market_vector) + length(truck_market_vector_D)), min = 0.00025, mode = 0.00049, max = 0.00075, shape = 4)
        
        beta_truck_loading_dist <- rpert((length(truck_loading_vector) + length(truck_loading_vector_D)), min = 0.0001, mode = 0.00027, max = 0.0004, shape = 4)
        
        beta_local_dist <- rpert((length(spatial_vector) + length(spatial_vector_D)), min = 0.0005, mode = 0.00075, max = 0.001, shape = 4)
        
      }
      
      if(length(number_accepted) == 10){
        
        beta_net_mode = ((6* mean(beta_accepted$net_mean)) - mean(beta_accepted$net_min) - mean(beta_accepted$net_max))/4
        
        beta_net_dist <- rpert((length(Pvector) + length(Pvector_D)), min = mean(beta_accepted$net_min), mode = beta_net_mode, max = mean(beta_accepted$net_max), shape = 4)
        
        beta_net_E_mode = ((6* mean(beta_accepted$net_E_mean)) - mean(beta_accepted$net_E_min) - mean(beta_accepted$net_E_max))/4
        
        beta_net_dist_E <- rpert(length(Pvector_E), min = mean(beta_accepted$net_E_min), mode = beta_net_E_mode, max = mean(beta_accepted$net_E_max), shape = 4)
        
        beta_truck_feed_mode = ((6* mean(beta_accepted$feed_truck_mean)) - mean(beta_accepted$feed_truck_min) - mean(beta_accepted$feed_truck_max))/4
        
        beta_truck_feed_dist <- rpert((length(truck_feed_vector) + length(truck_feed_vector_D)), min = mean(beta_accepted$feed_truck_min), mode = beta_truck_feed_mode, max = mean(beta_accepted$feed_truck_max), shape = 4)
        
        beta_truck_pig_mode = ((6* mean(beta_accepted$pig_truck_mean)) - mean(beta_accepted$pig_truck_min) - mean(beta_accepted$pig_truck_max))/4
        
        beta_truck_pig_dist <- rpert((length(truck_pig_vector) + length(truck_pig_vector_D)), min = mean(beta_accepted$pig_truck_min), mode = beta_truck_pig_mode, max = mean(beta_accepted$pig_truck_max), shape = 4)
        
        beta_truck_market_mode = ((6* mean(beta_accepted$market_truck_mean)) - mean(beta_accepted$market_truck_min) - mean(beta_accepted$market_truck_max))/4
        
        beta_truck_market_dist <- rpert((length(truck_market_vector) + length(truck_market_vector_D)), min = mean(beta_accepted$market_truck_min), mode = beta_truck_market_mode, max = mean(beta_accepted$market_truck_max), shape = 4)
        
        beta_truck_loading_mode = ((6* mean(beta_accepted$loading_truck_mean)) - mean(beta_accepted$loading_truck_min) - mean(beta_accepted$loading_truck_max))/4
        
        beta_truck_loading_dist <- rpert((length(truck_loading_vector) + length(truck_loading_vector_D)), min = mean(beta_accepted$loading_truck_min), mode = beta_truck_loading_mode, max = mean(beta_accepted$loading_truck_max), shape = 4)
        
        beta_local_mode = ((6* mean(beta_accepted$local_spread_mean)) - mean(beta_accepted$local_spread_min) - mean(beta_accepted$local_spread_max))/4
        
        beta_local_dist <- rpert((length(spatial_vector) + length(spatial_vector_D)), min = mean(beta_accepted$local_spread_min), mode = beta_local_mode, max = mean(beta_accepted$local_spread_max), shape = 4) 
        
      }
      
      if(length(number_accepted > 10)){
        
        k <- length(number_accepted)
        
        beta_net_mode = ((6* beta_accepted$net_mean[k]) - beta_accepted$net_min[k] - beta_accepted$net_max[k])/4
        
        beta_net_dist <- rpert((length(Pvector) + length(Pvector_D)), min = beta_accepted$net_min[k], mode = beta_net_mode, max = beta_accepted$net_max[k], shape = 4)
        
        beta_net_E_mode = ((6* beta_accepted$net_E_mean[k]) - beta_accepted$net_E_min[k] - beta_accepted$net_E_max[k])/4
        
        beta_net_dist_E <- rpert(length(Pvector_E), min = beta_accepted$net_E_min[k], mode = beta_net_E_mode, max = beta_accepted$net_E_max[k], shape = 4)
        
        beta_truck_feed_mode = ((6* beta_accepted$feed_truck_mean[k]) - beta_accepted$feed_truck_min[k] - beta_accepted$feed_truck_max[k])/4
        
        beta_truck_feed_dist <- rpert((length(truck_feed_vector) + length(truck_feed_vector_D)), min = beta_accepted$feed_truck_min[k], mode = beta_truck_feed_mode, max = beta_accepted$feed_truck_max[k], shape = 4)
        
        beta_truck_pig_mode = ((6* beta_accepted$pig_truck_mean[k]) - beta_accepted$pig_truck_min[k] - beta_accepted$pig_truck_max[k])/4
        
        beta_truck_pig_dist <- rpert((length(truck_pig_vector) + length(truck_pig_vector_D)), min = beta_accepted$pig_truck_min[k], mode = beta_truck_pig_mode, max = beta_accepted$pig_truck_max[k], shape = 4)
        
        beta_truck_market_mode = ((6* beta_accepted$market_truck_mean[k]) - beta_accepted$market_truck_min[k] - beta_accepted$market_truck_max[k])/4
        
        beta_truck_market_dist <- rpert((length(truck_market_vector) + length(truck_market_vector_D)), min = beta_accepted$market_truck_min[k], mode = beta_truck_market_mode, max = beta_accepted$market_truck_max[k], shape = 4)
        
        beta_truck_loading_mode = ((6* beta_accepted$loading_truck_mean[k]) - beta_accepted$loading_truck_min[k] - beta_accepted$loading_truck_max[k])/4
        
        beta_truck_loading_dist <- rpert((length(truck_loading_vector) + length(truck_loading_vector_D)), min = beta_accepted$loading_truck_min[k], mode = beta_truck_loading_mode, max = beta_accepted$loading_truck_max[k], shape = 4)
        
        beta_local_mode = ((6* beta_accepted$local_spread_mean[k]) - beta_accepted$local_spread_min[k] - beta_accepted$local_spread_max[k])/4
        
        beta_local_dist <- rpert((length(spatial_vector) + length(spatial_vector_D)), min = beta_accepted$local_spread_min[k], mode = beta_local_mode, max = beta_accepted$local_spread_max[k], shape = 4) 
        
      }
      
      #Transmission rates
      rateTransmission_net <- beta_net_dist[1:length(Pvector)] * Pvector #pig movements
      rateTransmission_net_D <- beta_net_dist[(length(Pvector) + 1):length(beta_net_dist)] * Pvector_D
      
      rateTransmission_truck_feed <- beta_truck_feed_dist[1:length(truck_feed_vector)] * truck_feed_vector #feed truck
      rateTransmission_truck_feed_D <- beta_truck_feed_dist[(length(truck_feed_vector) + 1):length(beta_truck_feed_dist)] * truck_feed_vector_D
      
      rateTransmission_truck_pig <- beta_truck_pig_dist[1:length(truck_pig_vector)] * truck_pig_vector #pig truck
      rateTransmission_truck_pig_D <- beta_truck_pig_dist[(length(truck_pig_vector) + 1):length(beta_truck_pig_dist)] * truck_pig_vector_D
      
      rateTransmission_truck_market <- beta_truck_market_dist[1:length(truck_market_vector)] * truck_market_vector #market truck
      rateTransmission_truck_market_D <- beta_truck_market_dist[(length(truck_market_vector) + 1):length(beta_truck_market_dist)] * truck_market_vector_D
      
      rateTransmission_truck_loading <- beta_truck_loading_dist[1:length(truck_loading_vector)] * truck_loading_vector #crew truck
      rateTransmission_truck_loading_D <- beta_truck_loading_dist[(length(truck_loading_vector) + 1):length(beta_truck_loading_dist)] * truck_loading_vector_D
      
      #Exposure rate
      rateexposure_net <- beta_net_dist_E * Pvector_E #pig movements from exposed farms
      
      #Spatial transmission rate
      spatial_vector <- spatial_vector[1,]
      
      spatial_vector_D <- spatial_vector_D[1,]
      
      rateTransmission_sp <- beta_local_dist[1:length(spatial_vector)] * spatial_vector
      rateTransmission_sp_D <- beta_local_dist[(length(spatial_vector) + 1):length(beta_local_dist)] * spatial_vector_D
      
      #Amount of incoming transmission through each route at each time step for each farm 
      trans_routes_aux <- data.frame(id = 1:length(movement_matrix_timestep[1, ]),
                                     day = i,
                                     t = epi_t,
                                     
                                     net = 1 - exp(-rateTransmission_net),
                                     net_D = 1 - exp(-rateTransmission_net_D),
                                     
                                     net_E = 1 - exp(-rateexposure_net),
                                     
                                     sp = 1 - exp(-rateTransmission_sp),
                                     sp_D = 1 - exp(-rateTransmission_sp_D),
                                     
                                     truck_feed = 1 - exp(-rateTransmission_truck_feed),
                                     truck_feed_D = 1 - exp(-rateTransmission_truck_feed_D),
                                     
                                     truck_pig = 1 - exp(-rateTransmission_truck_pig),
                                     truck_pig_D = 1 - exp(-rateTransmission_truck_pig_D),
                                     
                                     truck_market = 1 - exp(-rateTransmission_truck_market),
                                     truck_market_D = 1 - exp(-rateTransmission_truck_market_D),
                                     
                                     truck_loading = 1 - exp(-rateTransmission_truck_loading),
                                     truck_loading_D = 1 - exp(-rateTransmission_truck_loading_D))
      
      #Sum transmission rate 
      rateTransmission_E <- rateexposure_net +
        rateTransmission_sp +
        rateTransmission_sp_D +
        rateTransmission_truck_feed +
        rateTransmission_truck_feed_D +
        rateTransmission_truck_pig +
        rateTransmission_truck_pig_D +
        rateTransmission_truck_market +
        rateTransmission_truck_market_D +
        rateTransmission_truck_loading +
        rateTransmission_truck_loading_D
      
      rateTransmission_I <- 1 - exp(-rateTransmission_net)
      
      rateTransmission_E <- 1 - exp(-rateTransmission_E)
      
      rateTransmission_D <- 1 - exp(-rateTransmission_net_D)
      
      #biosecurity indication - biosecurity object is a vector of values for each farm
      ratehardlyinfected <- biosecurity * as.numeric(!I)
      rateTransmission_E <- rateTransmission_E - (rateTransmission_E * ratehardlyinfected)
      rateTransmission_I <- rateTransmission_I - (rateTransmission_I * ratehardlyinfected)
      rateTransmission_D <- rateTransmission_D - (rateTransmission_D * ratehardlyinfected)
      
      
      #Random probability a farm get infected
      Prand_I <- runif(n = length(movement_matrix_timestep[1, ]), min = 0, max = 1)
      Prand_E <- runif(n = length(movement_matrix_timestep[1, ]), min = 0, max = 1)
      Prand_D <- runif(n = length(movement_matrix_timestep[1, ]), min = 0, max = 1)
      
      #If rate of transmission is higher than the random probability of becoming infected then that farm will become exposed?
      PExposure <- as.numeric(rateTransmission_E >= Prand_E)
      PInfection <- as.numeric(rateTransmission_I >= Prand_I)
      PDetected <- as.numeric(rateTransmission_D >= Prand_D)
      
      #Adding newly infected farms to the exposed compartment
      Enew = as.numeric(S & PExposure & !PInfection & !PDetected) #farms which are susceptible and met the exposure threshold
      Inew = as.numeric((S & PInfection & !PDetected)|(E & PInfection & !PDetected)) #farms which are susceptible or exposed and met the infection threshold
      Dnew = as.numeric((S & PDetected)|(E & PDetected)|(I & PDetected))
      #looking for farms which are in both the susceptible and Pinfection group
      
      if(length(which(Enew == 1)) > 0){
        
        for (m in which(Enew == 1)){
          
          trans_true <- c(Pvector_E[m], truck_feed_vector[m], truck_feed_vector_D[m], truck_pig_vector[m], truck_pig_vector_D[m], truck_market_vector[m], truck_market_vector_D[m], truck_loading_vector[m], truck_loading_vector_D[m], spatial_vector[m], spatial_vector_D[m])
          
          trans_true <- as.numeric(trans_true > 0)
          
          b_E = data.frame(net_E = (beta_net_dist_E[m] * trans_true[1])/trans_true[1], 
                           feed_truck = (beta_truck_feed_dist[m] * trans_true[2] + beta_truck_feed_dist[(length(truck_feed_vector) + m)] * trans_true[3])/(trans_true[2] + trans_true[3]), 
                           pig_truck = (beta_truck_pig_dist[m] * trans_true[4] + beta_truck_pig_dist[(length(truck_pig_vector) + m)] * trans_true[5])/(trans_true[4] + trans_true[5]), 
                           market_truck = (beta_truck_market_dist[m] * trans_true[6] + beta_truck_market_dist[(length(truck_market_vector) + m)] * trans_true[7])/(trans_true[6] + trans_true[7]), 
                           loading_truck = (beta_truck_loading_dist[m] * trans_true[8] + beta_truck_loading_dist[(length(truck_loading_vector) + m)] * trans_true[9])/(trans_true[8] + trans_true[9]), 
                           local_spread = (beta_local_dist[m] * trans_true[10] + beta_local_dist[(length(spatial_vector) + m)] * trans_true[11])/(trans_true[10] + trans_true[11]))
          
          beta_values1_E <- rbind(beta_values1_E, b_E)
          
        }}
      
      if(length(which(Inew == 1)) > 0){  
        
        for(n in which(Inew == 1)){
          
          net2[length(net2) + 1] <- (beta_net_dist[n] * as.numeric(Pvector[n] > 0))/as.numeric(Pvector[n] > 0)
          
        }}
      
      if(length(which(Dnew == 1)) > 0){
        
        for(p in which(Dnew == 1)){
          
          net_D2[length(net_D2) + 1] <- (beta_net_dist[length(Pvector) + p] * as.numeric(Pvector_D[p] > 0))/as.numeric(Pvector_D[p] > 0)
          
        }}
      
      b_I <- data.frame(net = c(net2, net_D2))
      
      beta_values1_I <- rbind(beta_values1_I, b_I)
      
      
      M_Sim_E[, epi_t] = Enew #They are then added to the exposed group - may need to add a time exposed here
      M_Sim_I[, epi_t] = Inew
      
      M_Sim_D[, epi_t] = Dnew
      
      trans_routes_aux$exposed <- Enew #They are also added to the transmission route dataframe (why?)
      trans_routes_aux$infected <- Inew
      trans_routes_aux$detected <- Dnew
      trans_routes <- rbind(trans_routes, trans_routes_aux)
      
      #Randomly simulate the time the new farms will spend in the exposed compartment
      aux_tE = rep(0, length(Enew))
      if(length(aux_tE[Enew > 0]) > 0){
        aux_tE[Enew > 0] <- VGAM::rpospois(length(which(Enew != 0)), tE_time) 
      }
      
      #Moving farms into infected state once they have one day left in the exposed state
      M_Sim_I[which(tE == 1), epi_t] = 1
      M_Sim_I[which(M_Sim_I[, epi_t - 1] > 0), epi_t] = 1 +  M_Sim_I[which(M_Sim_I[, epi_t - 1] > 0), epi_t - 1]
      #M_Sim_I[which(M_Sim_I[, epi_t] > max_I), epi_t] <- 0
      
      M_Sim_I[which(I & PDetected), epi_t] <- 0
      M_Sim_I[which(detected_farm == 1), epi_t] = 0
      
      #So tI indicates how long they have until they leave the compartment? So this is reducing that time to indicate the passing of days?
      temp_tE = tE - 1
      tEnew = as.numeric(temp_tE < 0) + temp_tE
      
      #combining the current farms in these compartments (with time reduced by 1) and new farms to these compartments
      tE = tEnew + aux_tE
      
      #Add the new results to the infected and recovery matrix time
      M_Sim_E[, epi_t] = tE
      M_Sim_E[which((E & PInfection)|(E & PDetected)), epi_t] <- 0
      
      M_Sim_D[which(M_Sim_D[, epi_t - 1] > 0), epi_t] = 1 +  M_Sim_D[which(M_Sim_D[, epi_t - 1] > 0), epi_t - 1]
      M_Sim_D[which(detected_farm == 1), epi_t] <- 1
      #M_Sim_D[which(M_Sim_I[, epi_t - 1] == max_I), epi_t] <- 1
      
    }
    
    E_farm <- unique(which(M_Sim_E > 0, arr.ind = TRUE)[,1]) #farms which were exposed during the epidemic
    
    exp_row_col <- as.data.frame(which(M_Sim_E > 0, arr.ind = TRUE))
    
    E_data <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
    names(E_data) <- c("node_id", "t_step", "duration_E", "Sim")
    
    for (i in 1:length(E_farm)) {
      
      r <- which(exp_row_col$row == E_farm[i])[1]
      c <- exp_row_col$col[r]
      val <- M_Sim_E[E_farm[i], c]
      
      
      if(i < 2){
        
        E_data$node_id[i] <- E_farm[i]
        E_data$t_step[i] <- c
        E_data$duration_E[i] <- val
        E_data$Sim[i] <- w
        
      }else{
        E_vec <- c(E_farm[i], c, val, w)
        E_data <- rbind(E_data, E_vec)
      }
    }
    
    #Infected
    I_farm <- unique(which(M_Sim_I > 0, arr.ind = TRUE)[,1]) #farms which were infected during the epidemic
    
    inf_row_col <- as.data.frame(which(M_Sim_I > 0, arr.ind = TRUE))
    
    I_data <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
    names(I_data) <- c("node_id", "t_step", "duration_I", "Sim")
    
    for (i in 1:length(I_farm)) {
      
      r <- which(inf_row_col$row == I_farm[i])[1]
      c <- inf_row_col$col[r]
      val <- M_Sim_I[I_farm[i], c]
      
      if(i < 2){
        
        I_data$node_id[i] <- I_farm[i]
        I_data$t_step[i] <- c
        I_data$duration_I[i] <- val
        I_data$Sim[i] <- w
        
      }else{
        I_vec <- c(I_farm[i], c, val, w)
        I_data <- rbind(I_data, I_vec)
      }
    }
    
    #Detected
    D_farm <- unique(which(M_Sim_D > 0, arr.ind = TRUE)[,1]) #farms which became removed during the epidemic
    
    rem_row_col <- as.data.frame(which(M_Sim_D > 0, arr.ind = TRUE))
    
    D_data <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
    names(D_data) <- c("node_id", "t_step", "duration_D", "Sim")
    
    for (i in 1:length(D_farm)) {
      
      r <- which(rem_row_col$row == D_farm[i])[1]
      c <- rem_row_col$col[r]
      val <- M_Sim_D[D_farm[i], c]
      
      if(i < 2){
        
        D_data$node_id[i] <- D_farm[i]
        D_data$t_step[i] <- c
        D_data$duration_D[i] <- val
        D_data$Sim[i] <- w
        
      }else{
        D_vec <- c(D_farm[i], c, val, w)
        D_data <- rbind(D_data, D_vec)
      }
    }
    
    #transmission routes for infected farms
    trans_mod <- trans_routes %>% 
      filter(exposed == 1 | infected == 1) %>% 
      mutate(Sim = w)
    
    ### outbreaks
    sim_res_inf <- D_data %>% 
      mutate(farm_type = id_nodes$farm_type[node_id]) %>% 
      count(t_step, farm_type)
    
    day_results_inf_sow <- sim_res_inf[which(sim_res_inf$farm_type == "sow"),]
    day_results_inf_nur <- sim_res_inf[which(sim_res_inf$farm_type == "nursery"),]
    day_results_inf_fin <- sim_res_inf[which(sim_res_inf$farm_type == "finisher"),]
    
    sims_coords <- D_data %>% 
      mutate(lat = id_nodes$lat[node_id], long = id_nodes$long[node_id])
    
    dmx_sim <- distance_matrix(sims_coords$long,
                               sims_coords$lat,
                               cutoff = 900000,
                               min_dist = 100) #in meters
    
    # sow 
    simulated <- c(
      mu_total_sim_sow = sum(day_results_inf_sow$n),
      mu_median_sim_sow = mean(day_results_inf_sow$n),
      mu_max_sim_sow = max(day_results_inf_sow$n),
      
      mu_total_sim_fin = sum(day_results_inf_fin$n),
      mu_median_sim_fin = mean(day_results_inf_fin$n),
      mu_max_sim_fin = max(day_results_inf_fin$n),
      
      mu_total_sim_nur = sum(day_results_inf_nur$n),
      mu_median_sim_nur = mean(day_results_inf_nur$n),
      mu_max_sim_nur = max(day_results_inf_nur$n),
      
      mu_mean_dist_sim = mean(dmx_sim))
    
    simulated[which(is.na(simulated) == TRUE | simulated == "Inf" | simulated == "-Inf")] <- 0
    
    res2 <- c(w, simulated)
    
    accepted <- res2[2] < (observed[1] + error[1]) &
      res2[2] > (observed[1] - error[1]) &
      res2[3] < (observed[2] + error[2]) &
      res2[3] > (observed[2] - error[2]) &
      res2[4] < (observed[3] + error[3]) &
      res2[4] > (observed[3] - error[3]) &
      res2[5] < (observed[4] + error[4]) &
      res2[5] > (observed[4] - error[4]) &
      res2[6] < (observed[5] + error[5]) &
      res2[6] > (observed[5] - error[5]) &
      res2[7] < (observed[6] + error[6]) &
      res2[7] > (observed[6] - error[6]) &
      res2[8] < (observed[7] + error[7]) &
      res2[8] > (observed[7] - error[7]) &
      res2[9] < (observed[8] + error[8]) &
      res2[9] > (observed[8] - error[8]) &
      res2[10] < (observed[9] + error[9]) &
      res2[10] > (observed[9] - error[9]) &
      res2[11] < (observed[10] + error[10]) &
      res2[11] > (observed[10] - error[10])
    
    if(accepted == TRUE){
      
      if(length(number_accepted) == 0){
        
        number_accepted <- w
        
      }else{
        
        number_accepted <- append(number_accepted, w)
        
      }
      
      beta_summary <- data.frame(net_mean = mean(beta_values1_I$net, na.rm = TRUE),
                                 net_min = min(beta_values1_I$net, na.rm = TRUE),
                                 net_max = max(beta_values1_I$net, na.rm = TRUE),
                                 
                                 net_E_mean = mean(beta_values1_E$net_E, na.rm = TRUE),
                                 net_E_min = min(beta_values1_E$net_E, na.rm = TRUE),
                                 net_E_max = max(beta_values1_E$net_E, na.rm = TRUE),
                                 
                                 feed_truck_mean = mean(beta_values1_E$feed_truck, na.rm = TRUE),
                                 feed_truck_min = min(beta_values1_E$feed_truck, na.rm = TRUE),
                                 feed_truck_max = max(beta_values1_E$feed_truck, na.rm = TRUE),
                                 
                                 pig_truck_mean = mean(beta_values1_E$pig_truck, na.rm = TRUE),
                                 pig_truck_min = min(beta_values1_E$pig_truck, na.rm = TRUE),
                                 pig_truck_max = max(beta_values1_E$pig_truck, na.rm = TRUE),
                                 
                                 market_truck_mean = mean(beta_values1_E$market_truck, na.rm = TRUE),
                                 market_truck_min = min(beta_values1_E$market_truck, na.rm = TRUE),
                                 market_truck_max = max(beta_values1_E$market_truck, na.rm = TRUE),
                                 
                                 loading_truck_mean = mean(beta_values1_E$loading_truck, na.rm = TRUE),
                                 loading_truck_min = min(beta_values1_E$loading_truck, na.rm = TRUE),
                                 loading_truck_max = max(beta_values1_E$loading_truck, na.rm = TRUE),
                                 
                                 local_spread_mean = mean(beta_values1_E$local_spread, na.rm = TRUE),
                                 local_spread_min = min(beta_values1_E$local_spread, na.rm = TRUE),
                                 local_spread_max = max(beta_values1_E$local_spread, na.rm = TRUE))
      
      beta_accepted <- rbind(beta_accepted, beta_summary)
      
      if(length(number_accepted) == 1) {
        
        surveil_values_sow <- surveil_sow
        surveil_values_nur <- surveil_nur
        surveil_values_fin <- surveil_fin
        surveil_values_gilt <- surveil_gilt
        surveil_values_boar <- surveil_boar
        
      }
      
      if(length(number_accepted) > 1){
        
        surveil_values_sow <-  append(surveil_values_sow, surveil_sow)
        surveil_values_nur <- append(surveil_values_nur, surveil_nur)
        surveil_values_fin <- append(surveil_values_fin, surveil_fin)
        surveil_values_gilt <- append(surveil_values_gilt, surveil_gilt)
        surveil_values_boar <- append(surveil_values_boar, surveil_boar)
        
      }
      
      return(list(E_data = E_data, I_data = I_data, D_data = D_data, trans_mod = trans_mod, time_start = time_start, res2 = res2, beta_values_E = beta_values1_E, beta_values_I = beta_values1_I, accepted = 1, surveil_values = c(surveil_values_sow, surveil_values_nur, surveil_values_fin, surveil_values_gilt, surveil_values_boar), cutoff_val = sp_cutoff[dmx_num]))
      
    }else{
      
      return(list(E_data = E_data, I_data = I_data, D_data = D_data, trans_mod = trans_mod, time_start = time_start, res2 = res2, beta_values_E = beta_values1_E, beta_values_I = beta_values1_I, accepted = 0, surveil_values = c(surveil_sow, surveil_nur, surveil_fin, surveil_gilt, surveil_boar), cutoff_val = sp_cutoff[dmx_num]))
      
    }
}

close(pb)
stopCluster(cl)

save(ss, file = paste0("~/calibration/", ru, ".RData"))
}
