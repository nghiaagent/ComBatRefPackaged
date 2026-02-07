####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

####  Match quantiles
# keep_zero: zero values in the original counts don't change
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi, keep_zero=TRUE){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(keep_zero) {
        if(counts_sub[a, b] <= 1){
          new_counts_sub[a,b] <- counts_sub[a, b]
        }else{
          tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
          if(abs(tmp_p-1)<1e-4){
            new_counts_sub[a,b] <- counts_sub[a, b]  
            # for outlier count, if p==1, will return Inf values -> use original count instead
          }else{
            new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
          }
        }
      } else {
        tmp_p = pnbinom(counts_sub[a, b], mu=old_mu[a, b], size=1/old_phi[a])
  
        if(abs(tmp_p-1) < 1e-6) {
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        } else {
          new_counts_sub[a,b] <- qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}