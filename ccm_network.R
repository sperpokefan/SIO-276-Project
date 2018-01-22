library(rEDM)
library(parallel)
library(pbapply)
library(purrr)
library(dplyr)

# Compute CCM Network

compute_ccm_web <- function(in_file = "test/peb_clean.Rdata", best_E_file = "test/peb_best_E.csv",
                            num_surr = 500, num_cores = "automatic", 
                            surr_method = "random_shuffle", progress_bar = TRUE,
                            operating_system = .Platform$OS.type,
                            out_file = "test/peb_ccm_webs.Rdata", lib_auto = TRUE)
{
  #   for the ith row of ccm_params, do ccm calculations:
  #   use tp = 0 to measure CCM
  #   generate surrogate data
  #   compute significance of ccm vs surrogate data
  
  load(in_file)
  E_select <- read.csv(file = best_E_file, header = TRUE)
  
  # make sure best E is selected correctly
  
  taxa <- colnames(ts_matrix[,-1])
  E_select <- E_select[match(E_select$taxa,taxa),]
  
  drops_list <- which(is.na(E_select$E)) 
  drops <- E_select$taxa[drops_list] 
  
  ts_matrix <- ts_matrix[,!(names(ts_matrix) %in% drops)]
  E_select <- E_select[!(E_select$taxa %in% drops),]
  
  # variables to run all combinations
  n_obs <- NROW(ts_matrix)
  n_taxa <- NCOL(ts_matrix) 
  columns_to_use <- 2:n_taxa
  
  # generate all combinations where from != to
  ccm_params <- expand.grid(to = columns_to_use, from = columns_to_use)
  ccm_params <- ccm_params[ccm_params$from != ccm_params$to,]
  
  # code to run ccm for one pair
  run_ccm <- function(i)
  {
    temp_df <- ts_matrix[, c(ccm_params$from[i], ccm_params$to[i])]
    
    # select lib and pred
    # find first non-NA position: which(!is.na)
    # find max value of this position for from and to
    # start_value is after that: 1 +
    
    start_value <- 1 + max(min(which(!is.na(temp_df[, 1]))), 
                           min(which(!is.na(temp_df[, 2]))))
    
    if(lib_auto == TRUE){
      lib = c(start_value,n_obs)
      pred = lib
    }
    if(lib_auto == FALSE){
      lib = lib
      pred = pred
    }
    
    
    # Best E from s-map (hand-picked)
    
    best_E <- E_select$E[ccm_params$to[i]-1]
    
    # run ccm with tp = 0 and best_E
    
    ccm_tp_zero <- ccm(temp_df, lib = lib, pred = pred, lib_sizes = n_obs, 
                       E = best_E, tp = 0, lib_column = 1, target_column = 2, 
                       random_libs = FALSE, num_samples = 1, replace = FALSE,
                       first_column_time = FALSE, silent = TRUE)
    
    # make surrogate time series
    
    ts_index <- start_value:n_obs
    
    surrogate_data <- make_surrogate_data(as.numeric(temp_df[ts_index,2]), 
                                          method = surr_method, 
                                          num_surr = num_surr)
    
    # run ccm with surrogate data
    
    ccm_surrogate <- map_df(1:NCOL(surrogate_data), function(j) {
      temp_df[ts_index, 2] <- surrogate_data[,j]
      ccm(temp_df, lib = lib, pred = pred, lib_sizes = n_obs, 
          E = best_E, tp = 0, lib_column = 1, target_column = 2, 
          random_libs = FALSE, num_samples = 1, replace = FALSE,
          first_column_time = FALSE, silent = TRUE)
    })
    
    
    return(data_frame(from = ccm_params$from[i], 
                      to = ccm_params$to[i], 
                      best_E = best_E, 
                      ccm_tp_zero = list(ccm_tp_zero), 
                      ccm_surrogate = list(ccm_surrogate)))
  }
  
  # run in parallel
  
  # how many cores to use
  if(is.character(num_cores) == TRUE){no_cores <- detectCores() - 1
  } else {no_cores <- num_cores}
  
  # initiate cluster
  if(operating_system == "windows"){
    cluster <- makeCluster(no_cores)}
  if(operating_system != "windows"){
    cluster <- makeCluster(no_cores, type = "FORK")}
  
  # export variables
  clusterExport(cluster, varlist = c("ts_matrix", "E_select", "ccm_params", 
                                     "num_surr", "surr_method", "n_obs", "lib", "pred"),
                envir = environment())
  # export packages and run code
  if(progress_bar == TRUE){
    clusterEvalQ(cluster, expr =  c(library(rEDM),library(purrr), library(dplyr), library(pbapply)))
    ccm_results <- pblapply(cl = cluster, X = 1:NROW(ccm_params), FUN = run_ccm)}
  if(progress_bar == FALSE){
    clusterEvalQ(cluster, expr =  c(library(rEDM),library(purrr), library(dplyr)))
    ccm_results <- parLapply(cl = cluster, X = 1:NROW(ccm_params), fun = run_ccm)}
  
  # stop cluster
  stopCluster(cl = cluster)
  
  # convert to matrix
  n_tibble = length(ccm_results)
  
  # Create empty matrices to store results
  ccm_matrix = array(0,c(n_taxa,n_taxa)) # create an empty array to store results
  colnames(ccm_matrix) = names(ts_matrix)
  rownames(ccm_matrix) = names(ts_matrix)
  
  ccm_significance = ccm_matrix # create an empty array to store results
  bestE_matrix = ccm_matrix
  
  # Extract data in a loop
  for (i in 1:n_tibble){
    from = ccm_results[[i]]$from
    to = ccm_results[[i]]$to
    actual_rho = ccm_results[[i]]$ccm_tp_zero[[1]]$rho
    ccm_matrix[from,to] = actual_rho
    bestE_matrix[from,to] = ccm_results[[i]]$best_E
    #  Calculate and store p-value
    rho_surr = ccm_results[[i]]$ccm_surrogate[[1]]$rho
    ccm_significance[from,to] = (sum(actual_rho < rho_surr) + 1) / (length(rho_surr) + 1)
    
  }
  
  ccm_matrix = ccm_matrix[-1,-1]
  ccm_significance = ccm_significance[-1,-1]
  bestE_matrix = bestE_matrix[-1,-1]
  
  save(ccm_matrix, ccm_significance, bestE_matrix,file = out_file)
  
}