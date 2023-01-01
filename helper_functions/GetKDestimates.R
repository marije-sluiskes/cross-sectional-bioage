
################################################################################
### Function to obtain biological age estimates using the Klemera-Doubal method.
# NB: lower-case letter comments in function correspond to the steps 
# given in the Klemera-Doubal paper.
# ------------------------------------------------------------------------------
## IN: 
# - df_vars: dataframe with (candidate) biomarkers
# - c: chronological age vector
## OUT:
# - matrix with biological age estimates b_e (excluding chronological age as marker) and b_ec (including chronological age as marker)
################################################################################

GetKDestimates <- function(df_vars, c){
  
  numvar <- ncol(df_vars)
  n <- nrow(df_vars)
  
  # f
  mat_res <- matrix(nrow = numvar, ncol = 4)
  colnames(mat_res) <- c("a", "b", "s", "r") # intercept, slope, residual standard deviation, correlation coefficient
  
  for (i in 1:numvar){
    mod1 <- lm(df_vars[,i] ~ c)
    mat_res[i,1:2] <- mod1$coefficients
    mat_res[i,3] <- sigma(mod1)
    mat_res[i,4] <- cor(df_vars[,i], c)
  }
  mat_res <- as.data.frame(mat_res)
  
  # g
  r_char = ( sum( mat_res$r^2 / sqrt(1 - mat_res$r^2) ) ) /  sum( abs(mat_res$r) / sqrt(1 - mat_res$r^2) ) # slightly different from KD-paper, abs(r_j)
  s_char = ( sqrt(1-r_char^2) / r_char ) * ( (max(c) - min(c)) / sqrt(12) ) 
  
  # h
  b_e <- vector(length = nrow(df_vars))
  
  # loop over persons
  for(p in 1:n){
    
    num <- vector(length = numvar)
    denom <- vector(length = numvar)
    
    # loop over markers
    for (j in 1:numvar){
      num[j] = ( df_vars[p,j] - mat_res$a[j] ) * mat_res$b[j] / mat_res$s[j]^2
      denom[j] = mat_res$b[j]^2 / mat_res$s[j]^2
    } 
    
    b_e[p] = sum(num) / sum(denom)
  }
  
  diff_c_be <- b_e - c
  
  # s_be_sq is estimate of s_B, eq. 37
  s_be_sq <- ( sum( (diff_c_be - mean(diff_c_be))^2 ) / n ) - ( (1 - r_char^2) / r_char^2 ) * ( (max(c) - min(c))^2 / (12*numvar) )
  
  # i 
  b_ec <- vector(length = n)
  
  # loop over persons
  for(p in 1:n){
    
    num <- vector(length = numvar)
    denom <- vector(length = numvar)
    
    # loop over biomarkers
    for (j in 1:numvar){
      num[j] = ( df_vars[p,j] - mat_res$a[j] ) * mat_res$b[j] / mat_res$s[j]^2
      denom[j] = mat_res$b[j]^2 / mat_res$s[j]^2 
    } 
    
    b_ec[p] = ( sum(num) + c[p]/s_be_sq ) / ( (sum(denom) + 1/s_be_sq) )
  }
  
  return(cbind(b_e, b_ec))
  
}
