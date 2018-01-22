library(rEDM)
library(parallel)
library(pbapply)
library(purrr)
library(dplyr)

# Simplex/S-Map Analysis (Best-E & Nonlinearity)

simplex_smap_bestE <- function(in_file = "test/peb_clean.Rdata",
                               out_file = "test/peb_simplex_smap.Rdata",
                               Best_E_File = "test/peb_best_E.csv", text = FALSE, lib_auto = TRUE)
{
  
  load(in_file)
  
  data_output <- list()
  
  ts_matrix <- ts_matrix[,-1]
  
  # Selecting best E
  
  for (i in 1:length(ts_matrix[1,])) {
    
    # The length of your timeseries
    
    n <- length(ts_matrix[,i]) 
    ts <- as.numeric(ts_matrix[,i])
    
    start_value <- 1 + min(which(!is.na(ts)))
                           
    
    # use leave one out cross-validation (i.e., the whole time series is both library and predictor)
    
    
    if(lib_auto == TRUE){
    lib = c(start_value,n)
    pred = lib
    }
    if(lib_auto == FALSE){
      lib = lib
      pred = pred
    }
    
    E_list <- 1:10 # List of embedding dimensions to try in univariate_SSR
    rho_vs_e_output <- simplex(ts, lib = lib, pred = pred, E = E_list, silent = T)
    rho_vs_e_output$E <- E_list # Naming the first column so we know we are looking at embedding dimensions
    
    par(mar = c(5,4,3,3))
    layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
    
    plot(ts, type = "l", ylab = "Abundance", xlab = "Time", main = colnames(ts_matrix)[i])
    
    plot(rho_vs_e_output$E, rho_vs_e_output$rho,
         type = "l", xlab = "E", ylab = "rho", ylim = c(0,1), col = "blue")
    abline(h=0)
    
    
    if(isTRUE(text)==FALSE){
      best_E <- round(as.numeric(locator(n = 1)))[1]
    }
    
    if(isTRUE(text)==TRUE){
      abline(v = c(2,4,6,8), col = "lightgrey")
      best_E <- as.numeric(readline("E = "))
    }
    
    
    
    e_select <- list(rho_vs_e_output, best_E, ts)
    
    data_output[[i]] <- e_select
    
    names(data_output[[i]]) <- c("rho_vs_e_output", "Best_E", "ts")
    
    step <- paste0("E = ",best_E,", ",i,"/",length(ts_matrix[1,]))
    print(step)
    
  }
  
  for (i in 1:length(ts_matrix[1,])) {
    
    if(is.na(data_output[[i]]$Best_E)){
      data_output[[i]]$rho_vs_tp_output <- NA
      data_output[[i]]$rho_vs_theta_output <- NA
      data_output[[i]]$nonlinearity <- 0
      
    }
    if(!is.na(data_output[[i]]$Best_E)){
      
      if(data_output[[i]]$Best_E == 1){
        data_output[[i]]$rho_vs_tp_output <- NA
        data_output[[i]]$rho_vs_theta_output <- NA
        nonlinear <- 0
        data_output[[i]]$nonlinearity <- nonlinear
        
      }
      if(data_output[[i]]$Best_E > 1){
        # The length of your timeseries
        
        n <- length(ts_matrix[,i]) 
        ts <- as.numeric(ts_matrix[,i])
        
        # use leave one out cross-validation (i.e., the whole time series is both library and predictor)
        
        lib = lib
        pred = pred
        
        # Run rho vs. tp
        # ## Examining tp which looks at predictability the further you move into the future
        # ## things should normally be less predictable the further into the future you go
        
        tp_list <- 1:10 # number of steps ahead we are predicting 
        
        rho_vs_tp_output <- simplex(data_output[[i]]$ts, lib = lib, pred = pred, E = data_output[[i]]$Best_E, tp = tp_list,silent = TRUE)
        
        # Run rho vs. theta
        # remember that theta is looking at the relative weighting of the neighbors
        
        theta_list <- c(0,0.001,0.003,0.01,0.03,0.1,0.3,0.5,1,1.5,2,2.5,3,4,5,6,7,8) 
        rho_vs_theta_output <- s_map(data_output[[i]]$ts, lib = lib, pred = pred, E = data_output[[i]]$Best_E, theta = theta_list,silent = T)
        
        layout(matrix(c(1,1,2,2,3),1,5, byrow = TRUE))
        
        par(mar = c(4.5,4,5,2))
        plot(rho_vs_theta_output$theta, rho_vs_theta_output$rho,
             type = "l", xlab = "theta", ylab = "rho", col = "red")
        title(paste0(colnames(ts_matrix)[i],", E = ",data_output[[i]]$Best_E))
        
        par(mar = c(4.5,4,5,2))
        plot(rho_vs_tp_output$tp, rho_vs_tp_output$rho,
             type = "l", xlab = "tp", ylab = "rho", col = "green", ylim = c(0,1))
        
        par(mar = c(10,1,10,1))
        plot(1,1, xaxt = 'n', yaxt = 'n', pch = "Y", col = "green",
             xlim = c(-2,2), xlab = NA, ylab = NA, main = "Nonlinear?")
        points(-1,1, pch = "N", col = "red")
        abline(v = 0)
        
        
        if(isTRUE(text)==FALSE){
          nonlinear <- as.numeric(locator(n = 1))[[1]]
          if(nonlinear > 0){nonlinear = 1}else{nonlinear = 0}
        }
        
        if(isTRUE(text)==TRUE){nonlinear <- readline("Nonlinear? (If Yes, 1, if No, 0)  ")}
        
        if(nonlinear == 1){
          step <- paste0("Nonlinear? ","Yes",", ",i,"/",length(ts_matrix[1,]))
          print(step)
        }
        if(nonlinear == 0){
          step <- paste0("Nonlinear? ","No",", ",i,"/",length(ts_matrix[1,]))
          print(step)
        }
        
        
        data_output[[i]]$rho_vs_tp_output <- rho_vs_tp_output
        data_output[[i]]$rho_vs_theta_output <- rho_vs_theta_output
        data_output[[i]]$nonlinearity <- nonlinear
        
      }
      
    }
    
  }
  
  E_select <- as.data.frame(matrix(NA,nrow = length(ts_matrix[1,]),2))
  colnames(E_select) <- c("E", "taxa")
  
  for(i in 1:length(data_output)){
    
    if(!is.na(data_output[[i]]$Best_E == 1)){
      E_select[i,1] <- data_output[[i]]$Best_E
    }
    
    E_select[i,2] <- colnames(ts_matrix)[i]
    
  }
  
  save(data_output, file = out_file)
  write.csv(E_select, file = Best_E_File)
  
}
