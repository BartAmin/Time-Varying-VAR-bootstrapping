tvVAR_bootstrap = function(data, nIter, n, p){
  # --- Extract VAR parameters ----
  var_model = mvar(data = data,
                   type = rep("g", p), # Gaussian for each variable
                   level = rep(1, p), 
                   lags = 1, # Only lag 1 model
                   scale = FALSE,
                   pbar = FALSE)
  
  # Store intercept VAR parameters
  intercept <- unlist(var_model$intercepts)
  
  # Store phi parameters after converting sign if appropriate
  phi <- var_model$wadj[,,1]
  index = which(var_model$signs == -1)
  phi[index] = phi[index] * -1
  
  #Extract how many observations were included in model
  N = sum(var_model$call$data_lagged$included)
  
  # --- Peform parametric bootstrapping ----
  sdev <- array(NA, dim = c(p + 1, p, nIter))
  counter = 1
  while (counter <= nIter) {
    
    # Simulate stationary VAR data
    syn_data <- simulateVAR(pars = phi, 
                            means = intercept, 
                            Nt = N, 
                            residuals = 1) 
    
    tvvar_mod <- tryCatch({
      tvvarGAM(as.matrix(syn_data), 
               pbar = FALSE,
               thresholding = TRUE)
    }, error = function(e) {
      cat("Model fit error: Skipping iteration\n")
      return(NULL)  
    })
    
    # Skip iteration if model fitting failed
    if(is.null(tvvar_mod)){
      next
    }
    
    # Compute standard deviation of synthetic data across iterations
    sdev[,,counter] <- apply(tvvar_mod$Results_GAM$Estimate, c(1, 2), sd)
    
    # Print progress
    cat((counter / nIter) * 100, "%", "\n")
    
    counter = counter + 1
  }
  
  # --- Fit Time-Varying VAR on original data ----
  tvvar_mod = tvvarGAM(data,
                       thresholding = TRUE)
  
  # Compute standard deviation of original data across iterations
  real_sdev <- apply(tvvar_mod$Results_GAM$Estimate, c(1, 2), sd)
  
  # --- Calculate p values ----
  
  # Initialize matrix for p-values
  p_values <- matrix(NA, nrow = nrow(real_sdev), ncol = ncol(real_sdev))
  
  # Compute p-values using explicit loops
  for (j in 1:dim(real_sdev)[2]) {   
    for (i in 1:dim(real_sdev)[1]) { 
      p_values[i, j] <- mean(real_sdev[i, j] <= sdev[i, j, ])
    }
  }
  return(p_values)
}