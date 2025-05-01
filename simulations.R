# -----------------------------------------------------------------------
# Title: Time-Varying VAR Bootstrap Simulation
# Author: Bart Amin
# Description: Evaluates error control and ROC performance of p-values 
#              obtained via tvVAR bootstrapping. Compares stationary vs 
#              non-stationary stochastic processes
# -----------------------------------------------------------------------
library(mlVAR)
library(mgm)
library(devtools) 
library(tvvarGAM)
library(pROC)
library(ggplot2)
library(yardstick)
library(patchwork)


p_val_null = list()
p_val_var = list() 
nIter = 1000 # Amount of simulations in one iteration
nSim = 200 # Amount of iterations
n = 102 # Sample size (first two observations are omitted)
p = 5 # Number of variables in VAR

for(q in 1:nSim){
  
  # =====================
  # Stationary Bootstrap
  # =====================
  
  #--- Simulate Stationary VAR data ----
  stationary_data = matrix(0, nrow = n, ncol = p)  
  
  # Generate intercepts
  intercepts = runif(p, 3, 4)
  
  # Fill phi matrix with random fixed parameters, diagonal > off-diagonal
  phi_mat = matrix(runif(p * p, 0.01, 0.1) * sample(c(1, -1), (p * p), replace = TRUE), p, p)
  diag(phi_mat) = runif(p, 0.2, 0.5)
  
  # Simulate data
  for(i in 2:n){
    stationary_data[i,] = intercepts + phi_mat %*% stationary_data[i-1,] + rnorm(p)
  }
  stationary_data = stationary_data[-c(1:2),]
  
  # Perform tvVAR bootstrapping procedure to simulate p-values
  p_val_null[[q]] = tvVAR_bootstrap(stationary_data, nIter, n, p)
  
  # ========================
  # Non-stationary Bootstrap
  # ========================
  
  # --- Simulate Time Varying VAR data ----
  timevarying_data = matrix(0, nrow = n, ncol = p)  
  
  # Fill intercepts with random (increasing/decreasing) linear sequences
  intercepts = array(0, dim = c(1, p, n))
  for(i in 1:p){
    start = runif(1, 2.5, 3) 
    end = runif(1, 4, 4.5)
    if (sample(c(TRUE, FALSE), 1, prob = c(0.5, 0.5))) {
      direction = c(start, end)  # Increasing sequence
    } else {
      direction = c(end, start)  # Decreasing sequence
    }
    intercepts[,i,] = seq(direction[1], direction[2], length.out = n) 
  }
  
  # Fill phi matrix with random (increasing/decreasing) linear sequences, diagonal > off-diagonal
  phi_mat = array(0, dim = c(p, p, n))
  for(i in 1:p){
    for(j in 1:p){
      if(i == j){ # Phi's for diagonal
        phi_mat[i,j,] = seq(runif(1, 0.1, 0.2), runif(1, 0.5, 0.6), length.out = n) 
      } else{ # Phi's for off-diagonal
        phi_mat[i,j,] = seq(runif(1, 0.01, 0.02), runif(1, 0.09, 0.1), length.out = n) * 
                        sample(c(1, -1), 1, prob = c(0.5, 0.5))
      }
    }
  }
  
  for(i in 2:n){
    timevarying_data[i, ] = intercepts[,,i] + (phi_mat[,,i] %*% timevarying_data[i-1, ]) + rnorm(p)
  }
  timevarying_data = timevarying_data[-c(1:2),]
  
  # Perform tvVAR bootstrapping procedure to simulate p-values
  p_val_var[[q]] = tvVAR_bootstrap(timevarying_data, nIter, n, p)
}

# ==========================
# Evaluate type-1 error rate
# ==========================

# Set to your own path
#pdf("Error_control.pdf", width = 7, height = 5)

# ---- ECDF stationary p-values; evaluate type 1 error rate for each parameter ----

# Extract diagonal elements, excluding the first row (consisting of intercepts)
diagonal_values_stat <- sort(unlist(lapply(p_val_null, function(mat) diag(mat[-1, ]))))

# Extract intercepts
intercepts_stat <- sort(unlist(lapply(p_val_null, function(mat) mat[1, ])))

# Extract off-diagonal values
off_diagonal_values_stat <- sort(unlist(lapply(p_val_null, \(mat) {
trimmed <- mat[-1, ]           # Remove first row of intercepts
trimmed[!diag(nrow(trimmed))]  # Inlcude all off-diagonal elements
})))

df_stat <- data.frame(
Value = c(diagonal_values_stat, off_diagonal_values_stat, intercepts_stat),
Type = c(rep("Autoregressive", length(diagonal_values_stat)),
         rep("Cross-lagged", length(off_diagonal_values_stat)),
         rep("Intercepts", length(intercepts_stat)))
)

plot1 = ggplot(df_stat, aes(x = Value, color = Type)) +
          stat_ecdf(geom = "step") +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
          labs(title = paste('T = ', n - 2),   
               x = "P-values",
               y = "ECDF",
               color = "Type") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5)  
          )

# ============================
# Note on type-1 error rate
# ============================
# Only the p-values for auto-regressive effects have an ECDF that follows a uniform 
# distribution under the null hypothesis. Therefore we will only consider these effects 
# as only this type-1 error rate that is proportionate to a specified alpha level

# =================================================
# Calculate Receiver Operating Characteristic (ROC)
# =================================================

# Extract non-stationary p-values
diagonal_values_nonstat <- sort(unlist(lapply(p_val_var, function(mat) diag(mat[-1, ]))))

# Prepare data for ROC yardstick()
true = c(rep(1, length(diagonal_values_stat)), rep(0, length(diagonal_values_stat)))
prob = c(diagonal_values_nonstat, diagonal_values_stat)

df_roc <- data.frame(
  true = factor(true, levels = c(0, 1)),  
  prob  = prob
)

# Compute data for ROC curve
roc_data <- roc_curve(df_roc, true, prob)

plot2 = ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
          geom_line(color = 'red', size = 0.25) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
          labs(title = paste('T = ', n - 2),   
               x = "False Positive Rate",
               y = "True Positive Rate"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5)  
          )

# Print Area Under Curve (AUC)
print(roc_auc(df_roc, true, prob)$.estimate)

#Print both plots
plot1 + plot2

#dev.off()
