
################################################################################
### Function to get synthetic chronological age, biological age and markers.
# Data generated under the assumption that each marker X is governed by
# X = b_C * C + b_Delta * Delta + b_Lambda * Lambda + error
# ------------------------------------------------------------------------------
## IN: 
# - n: number of cases
# - c_min: minimum chronological age
# - c_max: maximum chronological age
# - s_delta: standard deviation of 'ageing rate' delta (and lambda)
# - beta_c: coefficient of chronological age
# - beta_delta: coefficient of delta 
# - beta_lambda: coefficient of lambda (= same distribution as delta, but not associated with biological age)
# - s_j: standard error of errors 
## OUT:
# - df_syn: matrix with length(beta_c) variables, chronological age c and biological age b
################################################################################

GenerateSyntheticData <- function(n, mean_c, sd_c, s_delta, beta_c, beta_delta, beta_lambda, s_j){

  m <- length(beta_c)
  
  c <- rnorm(n, mean_c, sd_c)
  delta <- rnorm(n, mean = 0, sd = s_delta)
  lambda <- rnorm(n, mean = 0, sd = s_delta)
  b <- c + delta
  
  x <- matrix(nrow = n, ncol = m)
  
  for (i in 1:m){
    x[,i]   <- beta_c[i] * c + beta_delta[i] * delta + beta_lambda[i] * lambda + rnorm(n, mean = 0, sd = s_j[i])
  }
  
  df_syn <- as.data.frame(cbind(x, c, b))
  
  return(df_syn)

}

