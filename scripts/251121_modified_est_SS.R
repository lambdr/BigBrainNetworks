# This script provides an alternative to the est_SS function in the MAPLE package. It is used to speed up the LDSC step 

# This is the core function that I modified 
est_SS <- function (exp_out = NULL, exposure = NULL, outcome = NULL, trait1.name = "exposure", trait2.name = "outcome", 
                     h2.fix.intercept = F, ldscore.dir = NULL, sumstats_dir = NULL, df_ldscore = NULL, M = 560836) {
  
### STEP 1: Data Preprocessing
  print("Loading and merging data")

  if (!is.null(exp_out)) {
    exp_out = stringr::str_split(exp_out, "_")[[1]]
    exposure = exp_out[1]
    outcome = exp_out[2]
  }
  
  # load data
  exp_raw <- load_raw_data(paste0(sumstats_dir, exposure, ".tsv.gz"))
  out_raw <- load_raw_data(paste0(sumstats_dir, outcome, ".tsv.gz"))
  
  # merge data with ldscores, filter out extreme values
  merged <- dplyr::left_join(
    df_ldscore, 
    exp_raw, 
    by = dplyr::join_by(SNP==variant_id)) %>% 
    dplyr::left_join(out_raw, by = dplyr::join_by(SNP==variant_id), suffix = c(".x", ".y")) %>% 
    dplyr::filter(abs(Z.x) < 5, abs(Z.y) < 5)
  
  merged <- merged %>% 
    dplyr::mutate(
      chi2.x = Z.x * Z.x,
      chi2.y = Z.y * Z.y,
      N.x = 4656, 
      N.y = N.x)
  
  
### STEP 2: Run ldsc_SS ###
# This section copies the necessary parts of MAPLE::ldsc_SS()
  print("Beginning data processing")
  Twostep=2 # I don't know why this works when set to 2 and not TRUE
  n.blocks=2

  cov <- matrix(NA, nrow = 2, ncol = 2)
  cov.se <- matrix(NA, nrow = 2, ncol = 2)
  I <- matrix(NA, nrow = 2, ncol = 2)
  I.se <- matrix(NA, nrow = 2, ncol = 2)
  merged <- merged[with(merged, order(CHR, BP)), ]
  n.snps = nrow(merged)
  merged$Zxy = merged$Z.x * merged$Z.y
  merged <- merged[, c("SNP", "chi2.x", "chi2.y", "N.x", "N.y", 
                       "Zxy", "L2")]
  merged11 = merged[, c("SNP", "chi2.x", "chi2.x", "N.x", 
                        "N.x", "chi2.x", "L2")]
  colnames(merged11) = c("SNP", "chi2.x", "chi2.y", "N.x", 
                         "N.y", "Zxy", "L2")
  merged22 = merged[, c("SNP", "chi2.y", "chi2.y", "N.y", 
                        "N.y", "chi2.y", "L2")]
  colnames(merged22) = c("SNP", "chi2.x", "chi2.y", "N.x", 
                         "N.y", "Zxy", "L2")
  
    step1.idx = which(merged$chi2.x < 30 & merged$chi2.y < 
                        30)
    seperator <- floor(seq(from = 1, to = length(step1.idx), 
                           length.out = (n.blocks + 1)))
    new.seperator <- c(1, step1.idx[seperator[2:n.blocks]], 
                       nrow(merged))

### STEP 3: Estimate heritability ### 
# This takes the most time of any step
    print("Beginning IRWLS estimates of h2")

    
    h2.1 = est_h2(merged = merged[, c("SNP", "chi2.x", "N.x", 
                                      "L2")], M, trait1.name, Twostep, h2.fix.intercept, 
                  step1.idx, n.blocks = n.blocks)
    h2.2 = est_h2(merged = merged[, c("SNP", "chi2.y", "N.y", 
                                      "L2")], M, trait2.name, Twostep, h2.fix.intercept, 
                  step1.idx, n.blocks = n.blocks)
    h2.1F = h2.1
    h2.2F = h2.2

### STEP 4: Est Gencov ### 
    print("Estimating gencov")
    h2_1 = est_gencov(merged = merged11, h1 = h2.1F$h2, 
                      h2 = h2.1F$h2, intercept.h1 = h2.1F$intercept, intercept.h2 = h2.1F$intercept, 
                      n.blocks = n.blocks, intercept.gencov = h2.1F$intercept, 
                      M, Twostep, fix.gcov.intercept = F, step1.idx, jknife = T)


    h2_2 = est_gencov(merged = merged22, h1 = h2.2F$h2, 
                      h2 = h2.2F$h2, intercept.h1 = h2.2F$intercept, intercept.h2 = h2.2F$intercept, 
                      n.blocks = n.blocks, intercept.gencov = h2.2F$intercept, 
                      M, Twostep, fix.gcov.intercept = F, step1.idx, jknife = T)

    rho_g = est_gencov(merged = merged, h1 = h2.1F$h2, h2 = h2.2F$h2, 
                       intercept.h1 = h2.1F$intercept, intercept.h2 = h2.2F$intercept, 
                       n.blocks = n.blocks, intercept.gencov = 0, M, Twostep, 
                       fix.gcov.intercept = F, step1.idx, jknife = T)

  I[1, 1] = h2_1$intercept
  I[2, 2] = h2_2$intercept
  I[1, 2] <- I[2, 1] <- rho_g$intercept
  
  ### END ###
  return(I)
}

### Supplemental Functions
# I wrote this function for data input
load_raw_data <- function(path){
  df <- data.table::fread(path) %>% 
    dplyr::mutate(Z = fixed_beta/fixed_sd) %>% 
    dplyr::select(variant_id, Z)
}

# These functions are copied directly from est_functions.R in MAPLE, with the only change
#   being that jackknife is disabled everywhere
##---------------------h2------------------------##
get.h2.weights <- function(h2, intercept, L2, N, M, x){
  h2 <- max(h2, 0)
  h2 <- min(h2, 1)
  wld <- as.numeric(lapply(X=L2, function(x){max(x,1)}))
  ld <- as.numeric(lapply(X=x, function(x){max(x,1)}))
  c <- h2*N/M
  het.w <- 1/(2*(intercept + c*ld )^2)
  oc.w <- 1/wld
  w <- sqrt(het.w*oc.w)
  return(w)
}


irwls_h2 <- function(y, L2, update.x, weights, intercept=1,
                     M, N, N.bar, fix.intercept=T, seperator){
  
  # calculate new weights
  weights = weights/sum(weights)
  x = L2 *N/N.bar
  
  if(fix.intercept){
    y = (y-intercept)
    for(i in 1:2){
      wy = y*weights
      wx = as.matrix(x*weights)
      fit = lm(wy~wx+0)
      h2 = drop(coef(fit)[1]) * M/N.bar
      h2 <- max(h2, 0)
      h2 <- min(h2, 1)
      # update weights
      weights = get.h2.weights(h2, intercept, L2, N, M, update.x)
      weights = weights/sum(weights)
    }
    ## preweight LD and chi:
    weighted.LD <- as.matrix(x*weights)
    weighted.chi <- as.matrix(y*weights)
  }else{
    for(i in 1:2){
      wy = y*weights
      wx = as.matrix(cbind(x*weights, weights))
      fit = lm(wy~wx+0)
      h2 = coef(fit)[1]*M/N.bar
      h2 <- max(h2, 0)
      h2 <- min(h2, 1)
      intercept = coef(fit)[2]
      weights = get.h2.weights(h2, intercept, L2, N, M, update.x)
      weights = weights/sum(weights)
    }
    ## preweight LD and chi:
    weighted.LD <- as.matrix(cbind(x*weights, weights))
    weighted.chi <- as.matrix(y*weights)
  }
  ## Perfrom analysis:
  n.blocks = length(seperator)-1
  n.snps <- length(y)
  p = ncol(weighted.LD)
  xty.block.values <- matrix(data=NA, nrow=n.blocks, ncol=p)
  xtx.block.values <- matrix(data=NA, nrow=p*n.blocks, ncol=p)
  
  from <- seperator
  to <- c(seperator[2:n.blocks]-1, n.snps)
  
  rep.from <- seq(from=1,to=p*n.blocks , by=p)
  rep.to <- seq(from =p,to=p*n.blocks ,by =p)
  
  colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
  
  for(i in 1:n.blocks){
    xty.block.values[i,] <- t(t(weighted.LD[from[i]:to[i],]) %*% weighted.chi[from[i]:to[i],])
    xtx.block.values[rep.from[i]:rep.to[i],] <- as.matrix(t(weighted.LD[from[i]:to[i],]) %*%
                                                            weighted.LD[from[i]:to[i],])
  }
  
  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(data=NA,nrow =p, ncol =p)
  colnames(xtx)<- colnames(weighted.LD)
  
  for(i in 1:p){
    xtx[i,] <- t(colSums(as.matrix(xtx.block.values[seq(from=i, to=p*n.blocks, by = p),])))
  }
  
  reg <- solve(xtx)%*% xty
  
  if(p==1){
    intercept.est <- intercept
    h2.est <- reg[1]/N.bar*M
  }
  
  if(p==2){
    intercept.est <- reg[2]
    h2.est <- reg[1]/N.bar*M
  }
  return(list(h2 = h2.est,
              intercept = intercept.est))
}


est_h2 <- function(merged, M, trait.name = NULL, Twostep=F, fix.intercept=T, step1.idx=NULL, n.blocks = 200){
  
  ## Read in summary statistics
  require(readr)
  if(is.null(trait.name)) trait.name = Trait
  colnames(merged) = c("SNP","chi2","N","L2")
  chi2 = merged$chi2
  N = merged$N
  L2 = merged$L2
  n.snps = nrow(merged)
  L2 <- as.numeric(lapply(X=L2, function(x){max(x,1)}))
  
  ## initial WEIGHTS:
  intercept = 1
  tot.agg <- (M*(mean(chi2)-1))/mean(L2*N)
  weights = get.h2.weights(tot.agg, intercept = 1, L2, N, M,L2)
  N.bar <- mean(N)
  x = L2 *N/N.bar
  
  if(fix.intercept) Twostep = F
  if(is.null(step1.idx)) Twostep = F
  if(Twostep){
    ### Two dtep estimator
    ## step I
    chi2.s1 <- chi2[step1.idx]
    L2.s1   <- L2[step1.idx]
    weights.s1 <- weights[step1.idx]
    N.s1 <- N[step1.idx]
    x.s1 <- x[step1.idx]
    seperator <- floor(seq(from=1, to=length(step1.idx), length.out =(n.blocks+1)))
    step1 = irwls_h2(chi2.s1, L2 = L2.s1, update.x=x.s1, weights.s1, intercept=1, M, N.s1, N.bar, fix.intercept, seperator)
    
    ## step 2
    new.seperator <-  c(1, step1.idx[seperator[2:n.blocks]], n.snps)
    step2 = irwls_h2(chi2,  L2=L2, update.x=L2, weights, intercept=step1$intercept, M, N, N.bar, fix.intercept=T, new.seperator)
    h2 = step2$h2
    intercept = step1$intercept
  }else{
    seperator <- floor(seq(from=1, to=length(chi2), length.out =(n.blocks+1)))
    step1 = irwls_h2(chi2,  L2=L2, update.x=L2, weights, intercept=1, M, N, N.bar, fix.intercept, seperator)
    h2 = step1$h2
    intercept = step1$intercept
  }
  return(list(h2 = h2,
              intercept = intercept))
}

##-------------------------gencov--------------------------##
get.gencov.weights <- function(L2, h1, h2, intercept.h1,
                               intercept.h2, intercept.gencov,
                               N.x, N.y, rho_g, N.bar, M, x){
  #### MAKE WEIGHTS:
  h1 <- max(h1, 0)
  h1 <- min(h1, 1)
  
  h2 <- max(h2, 0)
  h2 <- min(h2, 1)
  
  rho_g <- max(rho_g, -1)
  rho_g <- min(rho_g, 1)
  wld <- as.numeric(lapply(X=L2, function(x){max(x,1)}))
  ld <- as.numeric(lapply(X=x,function(x){max(x,1)}))
  a = N.x * h1 * ld/M + intercept.h1
  b = N.y * h2 * ld/M + intercept.h2
  c = sqrt(N.x*N.y) * rho_g * ld/M + intercept.gencov
  het.w <- 1/(a*b + c^2)
  oc.w <- 1/wld
  w <- sqrt(het.w * oc.w )
  return(w)
}


irwls_gencov <- function(y, L2, update.x, N,
                         weights,
                         N.x,N.y,
                         h1, h2,
                         intercept.h1,
                         intercept.h2,
                         intercept.gencov=0,
                         M, N.bar,
                         fix.intercept=T,
                         seperator) {
  
  jknife=F
  weights = weights/sum(weights)
  x = L2 *N/N.bar
  
  if(fix.intercept){
    y= y - intercept.gencov
    for(i in 1:2){
      wy = y*weights
      wx = as.matrix(x*weights)
      fit = lm(wy~wx+0)
      rho_g = drop(coef(fit)[1]* M/N.bar)
      rho_g <- max(rho_g, -1)
      rho_g <- min(rho_g, 1)
      # update weights
      weights = get.gencov.weights(L2, h1, h2,
                                   intercept.h1,
                                   intercept.h2,
                                   intercept.gencov,
                                   N.x, N.y,
                                   rho_g,
                                   N.bar, M, update.x)
      weights <-  weights/sum(weights)
    }
    ## preweight LD and chi:
    weighted.LD <- as.matrix(x*weights)
    weighted.chi <- as.matrix(y*weights)
  }else{
    for(i in 1:2){
      wy = y*weights
      wx = as.matrix(cbind(x*weights, 1*weights))
      fit = lm(wy~wx+0)
      rho_g = drop(coef(fit)[1]* M/N.bar)
      rho_g <- max(rho_g, -1)
      rho_g <- min(rho_g, 1)
      intercept.gencov = drop(coef(fit)[2])
      weights = get.gencov.weights(L2,h1, h2,
                                   intercept.h1,
                                   intercept.h2,
                                   intercept.gencov,
                                   N.x, N.y,
                                   rho_g,
                                   N.bar, M, update.x)
      weights <-  weights/sum(weights)
    }
    ## preweight LD and chi:
    weighted.LD <- as.matrix(cbind(x*weights, 1*weights))
    weighted.chi <- as.matrix(y*weights)
  }
  
  n.blocks = length(seperator)-1
  n.snps <- length(y)
  p = ncol(weighted.LD)
  xty.block.values <- matrix(data=NA, nrow=n.blocks, ncol=p)
  xtx.block.values <- matrix(data=NA, nrow=p*n.blocks, ncol=p)
  from <- seperator
  to <- c(seperator[2:n.blocks]-1, n.snps)
  rep.from <- seq(from=1,to=p*n.blocks , by=p)
  rep.to <- seq(from =p,to=p*n.blocks ,by =p)
  colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
  
  
  for(i in 1:n.blocks){
    xty.block.values[i,] <- t(t(weighted.LD[from[i]:to[i],]) %*% weighted.chi[from[i]:to[i],])
    xtx.block.values[rep.from[i]:rep.to[i],] <- as.matrix(t(weighted.LD[from[i]:to[i],]) %*%
                                                            weighted.LD[from[i]:to[i],])
  }
  
  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(data=NA,nrow =p, ncol =p)
  colnames(xtx)<- colnames(weighted.LD)
  
  for(i in 1:p){
    xtx[i,] <- t(colSums(as.matrix(xtx.block.values[seq(from=i, to=p*n.blocks, by = p),])))
  }
  
  reg <- solve(xtx)%*% xty
  
  if(jknife){
    
    # perform jacknife
    delete.from <- seq(from=1,to= p* n.blocks, by=p)
    delete.to <- seq(from=p,  to=p* n.blocks, by=p)
    delete.values <- matrix(data=NA, nrow=n.blocks, ncol = p)
    colnames(delete.values)<- colnames(weighted.LD)
    
    for(i in 1:n.blocks){
      xty.delete <- xty-xty.block.values[i,]
      xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
      delete.values[i,] <- solve(xtx.delete)%*% xty.delete
    }
    
    delete.values <- as.matrix(delete.values[,1:p])
    
    
    pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=p)
    colnames(pseudo.values)<- colnames(weighted.LD)
    for(i in 1:n.blocks){
      pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])
    }
    
    jackknife.cov <- cov(pseudo.values)/n.blocks
    jackknife.se <- sqrt(diag(jackknife.cov))
    
    coef.cov <- jackknife.cov[1,1]/(N.bar^2)*M^2
    rho_g.se <- sqrt(coef.cov)
    
    if(p==1){
      intercept.est <- intercept.gencov
      rho_g.est <- reg[1]/N.bar*M
      intercept.se <- NA
    }
    
    if(p==2){
      
      intercept.est <- reg[2]
      rho_g.est <- reg[1]/N.bar*M
      intercept.se <- jackknife.se[length(jackknife.se)]
      
    }
    return(list(rho_g = rho_g.est, rho_g.se= rho_g.se,
                intercept = intercept.est,
                intercept.se = intercept.se,
                delete.values =  delete.values,
                jk.est = reg[1]))
  }else{
    if(p==1){
      intercept.est <- intercept.gencov
      rho_g.est <- reg[1]/N.bar*M
    }
    
    if(p==2){
      intercept.est <- reg[2]
      rho_g.est <- reg[1]/N.bar*M
    }
    return(list(rho_g = rho_g.est,
                intercept = intercept.est))
  }
}


est_gencov <- function(merged,
                       h1, h2,
                       intercept.h1,
                       intercept.h2,
                       intercept.gencov=0,
                       M,Twostep=F,
                       fix.gcov.intercept=F,
                       step1.idx=NULL,
                       n.blocks=200,
                       jknife=F){
  
  Zxy = merged$Zxy
  L2 = merged$L2
  N.x = merged$N.x
  N.y = merged$N.y
  N = sqrt(N.x * N.y)
  N.bar <- mean(N)
  x = L2 *N/N.bar
  rho_g <- mean(Zxy)*M/mean(N*L2)
  
  weights = get.gencov.weights(L2, h1, h2,
                               intercept.h1,
                               intercept.h2,
                               intercept.gencov,
                               N.x, N.y,
                               rho_g, N.bar, M, L2)
  
  if(fix.gcov.intercept)  Twostep=F
  if(is.null(step1.idx)) Twostep = F
  if(Twostep==T){
    #message("Using two-step estimator with cutoff at 30.")
    ### Two step estimator
    Zxy.s1 <-Zxy[step1.idx]
    L2.s1   <- L2[step1.idx]
    weights.s1 <- weights[step1.idx]
    N.s1 <- N[step1.idx]
    x.s1 <- x[step1.idx]
    N.x.s1 <- N.x[step1.idx]
    N.y.s1 <- N.y[step1.idx]
    n.snps <- nrow(merged)
    seperator <- floor(seq(from=1, to=length(step1.idx), length.out =(n.blocks+1)))
    
    step1 = irwls_gencov( y = Zxy.s1,
                          L2 = L2.s1,
                          update.x = x.s1,
                          N=N.s1,
                          weights = weights.s1,
                          N.x = N.x.s1,
                          N.y = N.y.s1,
                          h1, h2,
                          intercept.h1,
                          intercept.h2,
                          intercept.gencov=0,
                          M, N.bar,
                          fix.intercept=F,
                          seperator)
    
    
    ## step 2
    #cat("step II \n")
    new.seperator <-  c(1, step1.idx[seperator[2:n.blocks]], n.snps)
    
    step2 =irwls_gencov(y=Zxy,
                        L2=L2,
                        update.x = L2,
                        N = N,
                        weights = weights,
                        N.x = N.x,
                        N.y = N.y,
                        h1, h2,
                        intercept.h1,
                        intercept.h2,
                        intercept.gencov=step1$intercept,
                        M, N.bar,
                        fix.intercept=T,
                        new.seperator)
    intercept = step1$intercept
    
    if(jknife){
      ## cpmbine step 1 and step 2
      x = L2 *N/N.bar
      c = sum(weights^2 *x)/sum(weights^2*x^2)
      est = c(step2$jk.est, step1$intercept)
      delete.values = matrix(0, nrow = n.blocks, ncol = 2)
      delete.values[,2] = step1$delete.values[,2]
      delete.values[,1] = M/N.bar *(step2$delete.values - c * (step1$delete.values[,2] - step1$intercept))
      
      pseudo.values = matrix(n.blocks*est,nrow =n.blocks,ncol=2, byrow = T) - (n.blocks-1)* delete.values
      
      jackknife.cov <- cov(pseudo.values)/n.blocks
      jackknife.se <- sqrt(diag(jackknife.cov))
      
      intercept.se <- jackknife.se[2]
      
    }
    
  }else{
    
    seperator <- floor(seq(from=1, to=nrow(merged), length.out =(n.blocks+1)))
    
    step1 = irwls_gencov( y = Zxy,
                          L2 = L2,
                          update.x = L2,
                          N=N,
                          weights = weights,
                          N.x = N.x,
                          N.y = N.y,
                          h1, h2,
                          intercept.h1,
                          intercept.h2,
                          intercept.gencov,
                          M, N.bar,
                          fix.intercept=fix.gcov.intercept,
                          seperator)
    intercept = step1$intercept
    if(jknife){
      intercept.se = step1$intercept.se
    }
    
  }
  
  if(jknife){
    return(list(intercept = intercept,
                intercept.se = intercept.se))
  }else{
    return(list(intercept = intercept))
  }
  
  
}
