library(rjags)
library(doParallel)
library(foreach)
#library(peakRAM)

dc_melding_merging_jags <- function(data,
                                    alpha0,
                                    alpha2,
                                    rho,
                                    Ntotal,
                                    burn_in) {
  tt <- length(data)
  
  model <- "model {
    for (i in 1:50) {
    p[i] <- 1/50
    }
    
    logit(deltaJF) <- alpha0
    logit(deltaAF) <- alpha0 + alpha2
    
    xJ[1] ~ dcat(p[1:50])
    sur[1] ~ dcat(p[1:50])
    imm[1] ~ dcat(p[1:50])
    
    alpha0 ~ dnorm(0, 1/4) T(-10, 10)
    alpha2 ~ dnorm(0, 1/4) T(-10, 10)
    alpha6 ~ dnorm(0, 1/4) T(-10, 10)
    eta_t <- exp(alpha6)
    
    rho ~ dunif(0, 10)
    
    for (tt in 2:t) {
      rateJ[tt-1] <- .5 * rho * deltaJF * x[tt-1]
      xJ[tt] ~ dpois(rateJ[tt-1])
      sur[tt] ~ dbin(deltaAF, x[tt-1])
      rate_imm[tt-1] <- x[tt-1] * eta_t
      imm[tt] ~ dpois(rate_imm[tt-1])
    }
    
    for (tt in 1:t) {
      x[tt] <- xJ[tt] + sur[tt] + imm[tt]
      y[tt] ~ dpois(x[tt])
    }
  }"
  
  foo <- jags.model(textConnection(model),
                    data = list(y = data,
                                t = tt,
                                alpha0 = alpha0,
                                alpha2 = alpha2,
                                rho = rho),
                    quiet = T)
  
  update(foo, burn_in, progress.bar = "none")
  
  jags <- coda.samples(model = foo,
                       variable.names = c("alpha6", "xJ", "sur", "imm"),
                       n.iter = Ntotal,
                       progress.bar = "none")
  out_jags <- as.matrix(jags)
  
  alpha6_jags <- out_jags[,1]
  imm_jags <- out_jags[,2:27]
  sur_jags <- out_jags[,28:53]
  xJ_jags <- out_jags[,54:79]
  
  return(list(alpha6_jags = alpha6_jags,
              imm_jags = imm_jags,
              sur_jags = sur_jags,
              xJ_jags = xJ_jags))
}


dc_melding_smc_jags <- function(data,
                                alpha0,
                                alpha2,
                                rho,
                                psi2_inits,
                                alpha_j,
                                N,
                                Ntotal,
                                n_chains,
                                toggle = FALSE) {
  tt <- length(data)
  
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  alpha0 <- matrix(alpha0, ncol = n_chains)
  alpha2 <- matrix(alpha2, ncol = n_chains)
  rho <- matrix(rho, ncol = n_chains)
  
  
  model <- "model {
  
    for (i in 1:50) {
      p[i] <- 1/50
    }
  
    for (i in 1:N) {
      logit(deltaJF[i]) <- alpha0[i]
      logit(deltaAF[i]) <- alpha0[i] + alpha2[i]  
  
      xJ[i,1] ~ dcat(p[1:50])
      sur[i,1] ~ dcat(p[1:50])
      imm[i,1] ~ dcat(p[1:50])
     
      # alpha0[i] ~ dnorm(0, 1/4) T(-10, 10)
      # alpha2[i] ~ dnorm(0, 1/4) T(-10, 10)
      alpha6[i] ~ dnorm(0, 1/4) T(-10, 10)
      eta_t[i] <- exp(alpha6[i])
       
      # rho[i] ~ dunif(0, 10)
       
      for (tt in 2:t) {
        rateJ[i,tt-1] <- .5 * rho[i] * deltaJF[i] * x[i,tt-1]
        xJ[i,tt] ~ dpois(rateJ[i,tt-1])
        sur[i,tt] ~ dbin(deltaAF[i], x[i,tt-1])
        rate_imm[i,tt-1] <- x[i,tt-1] * eta_t[i]
        imm[i,tt] ~ dpois(rate_imm[i,tt-1])
      }
    
      for (tt in 1:t) {
        x[i,tt] <- xJ[i,tt] + sur[i,tt] + imm[i,tt]
        logL[i,tt] <- logdensity.pois(y[i,tt], x[i,tt])
        phi[i,tt] <- -alpha * logL[i,tt] + 1e20
        zeros[i,tt] ~ dpois(phi[i,tt])
      }
    }
  }"
  
  
  alpha6_0 <- matrix(psi2_inits$alpha6_inits, ncol = n_chains)
  xJ_0 <- array(rep(10, tt*N), dim = c(tt, N/n_chains, n_chains))
  xJ_0 <- aperm(xJ_0, c(2,1,3))
  sur_0 <- array(rep(10, tt*N), dim = c(tt, N/n_chains, n_chains))
  sur_0 <- aperm(sur_0, c(2,1,3))
  imm_0 <- array(rep(10, tt*N), dim = c(tt, N/n_chains, n_chains))
  imm_0 <- aperm(imm_0, c(2,1,3))
  
  inits <- list(alpha6_inits = alpha6_0,
                xJ_inits = xJ_0,
                sur_inits = sur_0,
                imm_inits = imm_0)
  
  
  #tim_init <- Sys.time()
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "alpha0", "alpha2", "rho", "inits", "model", "toggle"),
                envir = environment())
  
  #mcmc_list <- foreach(i=1:n_chains, .packages = c("peakRAM","rjags")) %dopar% {}
  mcmc_list <- foreach(i=1:n_chains, .packages = c("rjags")) %dopar% {
    #peakRAM({
    if (toggle) {
      foo <- jags.model(textConnection(model),
                        data = list(y = DATA,
                                    t = tt,
                                    N = N/n_chains,
                                    alpha0 = alpha0[,i],
                                    alpha2 = alpha2[,i],
                                    rho = rho[,i],
                                    alpha = alpha_j,
                                    zeros = matrix(rep(0, tt*N/n_chains), ncol=tt)))
    } else {
      foo <- jags.model(textConnection(model),
                        data = list(y = DATA,
                                    t = tt,
                                    N = N/n_chains,
                                    alpha0 = alpha0[,i],
                                    alpha2 = alpha2[,i],
                                    rho = rho[,i],
                                    alpha = alpha_j,
                                    zeros = matrix(rep(0, tt*N/n_chains), ncol=tt)),
                        inits = list(alpha6 = inits$alpha6_inits[,i],
                                     xJ = inits$xJ_inits[,,i],
                                     sur = inits$sur_inits[,,i],
                                     imm = inits$imm_inits[,,i]))
    }
    
    
    jags <- coda.samples(model = foo,
                         variable.names = c("alpha6", "xJ", "sur", "imm"),
                         n.iter = Ntotal)
    #})
    
  }
  
  stopCluster(cl)
  
  # tim_end <- Sys.time()
  # tim_run <- tim_end - tim_init
  
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("alpha6", "xJ", "sur", "imm")
  var_names_index <- 1:26
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  var_names_full <- rep(NA, 79)
  var_names_full[1] <- paste0("^", var_names[1], "\\[")
  for (i in 1:tt) {
    var_names_full[i+1] <- paste0("^", var_names[2], "\\[[0-9]+,", var_names_index[i], "\\]")
  }
  for (i in 1:tt) {
    var_names_full[i+tt+1] <- paste0("^", var_names[3], "\\[[0-9]+,", var_names_index[i], "\\]")
  }
  for (i in 1:tt) {
    var_names_full[i+2*tt+1] <- paste0("^", var_names[4], "\\[[0-9]+,", var_names_index[i], "\\]")
  }
  
  alpha6_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  alpha6_res <- unlist(alpha6_res)
  alpha6_res <- matrix(alpha6_res, nrow = Ntotal)[Ntotal,]
  
  xJ_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    xJ_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+1], colnames(jags_out1))])
    xJ_temp <- unlist(xJ_temp)
    xJ_res[,i] <- matrix(xJ_temp, nrow = Ntotal)[Ntotal,]
  }
  
  sur_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    sur_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+tt+1], colnames(jags_out1))])
    sur_temp <- unlist(sur_temp)
    sur_res[,i] <- matrix(sur_temp, nrow = Ntotal)[Ntotal,]
  }
  
  imm_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    imm_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+2*tt+1], colnames(jags_out1))])
    imm_temp <- unlist(imm_temp)
    imm_res [,i] <- matrix(imm_temp, nrow = Ntotal)[Ntotal,]
  }
  
  
  return(list(alpha6_res = alpha6_res,
              xJ_res = xJ_res,
              sur_res = sur_res,
              imm_res = imm_res))
}









recap_jags <- function(data,
                       alpha_j,
                       inits,
                       N,
                       Ntotal,
                       n_chains) {
  data_f_j <- data$recap_f_j 
  data_f_a <- data$recap_f_a 
  data_m_j <- data$recap_m_j 
  data_m_a <- data$recap_m_a
  
  tt <- dim(data_f_j)[2]
  
  DATA_f_j <- array(rep(data_f_j, N/n_chains), 
                    dim = c(nrow(data_f_j), ncol(data_f_j), N/n_chains))
  DATA_f_a <- array(rep(data_f_a, N/n_chains), 
                    dim = c(nrow(data_f_a), ncol(data_f_a), N/n_chains))
  DATA_m_j <- array(rep(data_m_j, N/n_chains), 
                    dim = c(nrow(data_m_j), ncol(data_m_j), N/n_chains))
  DATA_m_a <- array(rep(data_m_a, N/n_chains), 
                    dim = c(nrow(data_m_a), ncol(data_m_a), N/n_chains))
  
  model <- "model{
    for (i in 1:N) {
      #----------Part 1. Regression equations (44)-------------#
      ##------------------------------------------------------------##
      ## alpha = (alpha0, alpha1, alpha2, alpha4)
      ##------------------------------------------------------------##
      for (t in 1:(tt-1)) {
        # Recapture rate pi
        logit(pM[i,t]) <- alpha[i,4] + alpha5[i,t]   # Male
        logit(pF[i,t]) <- alpha5[i,t]              # Female
      }
      
      # Juvenile survival rate deltaj
      logit(deltajM[i]) <- alpha[i,1] + alpha[i,2]  # Male
      logit(deltajF[i]) <- alpha[i,1]             # Female  
    
      # Adult survival rate deltaa
      logit(deltaaM[i]) <- alpha[i,1] + alpha[i,2] + alpha[i,3]   # Male
      logit(deltaaF[i]) <- alpha[i,1] + alpha[i,3]              # Female
      
      #----------Part 2. Priors on alpha's-------------#
      alpha[i,1] ~ dunif(-1e9, 1e9)
      alpha[i,3] ~ dunif(-1e9, 1e9)
    
      alpha[i,2] ~ dnorm(0, 1/4) T(-10, 10)
      alpha[i,4] ~ dnorm(0, 1/4) T(-10, 10)
    
      for (t in 1:(tt-1)) {
        alpha5[i,t] ~ dnorm(1, 1/4) T(-10, 10)
      }
      
      #----------Part 3. Likelihood-------------#
      ##------------------------------------------------------------##
      ## Female capture recapture rate
      ##------------------------------------------------------------##
      ### Rate for juveniles ###
      for (t in 1:(tt-1)) {
        logLjF[i,t] <- logdensity.multi(MjF[t, 1:tt, i], QjF[t, 1:tt, i], RjF[i,t])
        phijF[i,t] <- -alpha_j * logLjF[i,t] + 1000
        zerosjF[i,t] ~ dpois(phijF[i,t])
        qF[i,t] <- 1 - pF[i,t]  # 1 - pi
      }
      
      #### Definition on QjF ####
      for (t in 1:(tt-1)) {
        # Main diagonal
        QjF[t, t, i] <- deltajF[i] * pF[i,t]
        
        # Above main diagonal
        for (j in (t+1):(tt-1)) {
          QjF[t, j, i] <- deltajF[i] * deltaaF[i] * pF[i,j] * prod(qF[i, t:(j-1)]) # Is deltajF^j or just deltajF?
        }
      
        # Below main diagonal
        for (j in 1:(t-1)) {
          QjF[t, j, i] <- 0
        }
      
        # Last column
        QjF[t, tt, i] <- 1 - sum(QjF[t, 1:(tt-1), i])
      }
      
      ### Rate for adults ###
      for (t in 1:(tt-1)) {
        logLaF[i,t] <- logdensity.multi(MaF[t, 1:tt, i], QaF[t, 1:tt, i], RaF[i,t])
        phiaF[i,t] <- -alpha_j * logLaF[i,t] + 1000
        zerosaF[i,t] ~ dpois(phiaF[i,t])
      }
      
      #### Definition on QaF ####
      for (t in 1:(tt-1)) {
        # Main diagonal
        QaF[t, t, i] <- deltaaF[i] * pF[i,t]
      
        # Above main diagonal
        for (j in (t+1):(tt-1)) {
          QaF[t, j, i] <- deltaaF[i] * pF[i,j] * prod(qF[i, t:(j-1)]) # Is deltaaF^j or just deltaaF?
        }
      
        # Below main diagonal
        for (j in 1:(t-1)) {
          QaF[t, j, i] <- 0
        }
      
        # Last column
        QaF[t, tt, i] <- 1 - sum(QaF[t, 1:(tt-1), i])
      }
      
      
      ##------------------------------------------------------------##
      ## Male capture recapture rate
      ##------------------------------------------------------------##
      ### Rate for juveniles ###
      for (t in 1:(tt-1)) {
        logLjM[i,t] <- logdensity.multi(MjM[t, 1:tt, i], QjM[t, 1:tt, i], RjM[i,t])
        phijM[i,t] <- -alpha_j * logLjM[i,t] + 1000
        zerosjM[i,t] ~ dpois(phijM[i,t])
        qM[i,t] <- 1 - pM[i,t]
      }
      
      #### Definition on QjM ####
      for (t in 1:(tt-1)) {
        # Main diagonal
        QjM[t, t, i] <- deltajM[i] * pM[i,t]
        
        # Above main diagonal
        for (j in (t+1):(tt-1)) {
          QjM[t, j, i] <- deltajM[i] * deltaaM[i] * pM[i,j] * prod(qM[i, t:(j-1)]) # Is deltajM^j or just deltajF?
        }
        
        # Below main diagonal
        for (j in 1:(t-1)) {
          QjM[t, j, i] <- 0
        }
        
        # Last column
        QjM[t, tt, i] <- 1 - sum(QjM[t, 1:(tt-1), i])
      }
      
      ### Rate for adults ###
      for (t in 1:(tt-1)) {
        logLaM[i,t] <- logdensity.multi(MaM[t, 1:tt, i], QaM[t, 1:tt, i], RaM[i,t])
        phiaM[i,t] <- -alpha_j * logLaM[i,t] + 1000
        zerosaM[i,t] ~ dpois(phiaM[i,t])
      }
      
      #### Definition on QaF ####
      for (t in 1:(tt-1)) {
        # Main diagonal
        QaM[t, t, i] <- deltaaM * pM[i,t]
        
        # Above main diagonal
        for (j in (t+1):(tt-1)) {
          QaM[t, j, i] <- deltaaM[i] * pM[i,j] * prod(qM[i, t:(j-1)]) # Is deltaaM^j or just deltaaF?
        }
        
        # Below main diagonal
        for (j in 1:(t-1)) {
          QaM[t, j, i] <- 0
        }
        
        # Last column
        QaM[t, tt, i] <- 1 - sum(QaM[t, 1:(tt-1), i])
      }
    }
  }"
  
  
  alpha0_0 <- matrix(inits$alpha0_0, ncol = n_chains)
  alpha2_0 <- matrix(inits$alpha2_0, ncol = n_chains)
  alpha1_0 <- matrix(inits$alpha1_0, ncol = n_chains)
  alpha4_0 <- matrix(inits$alpha4_0, ncol = n_chains)
  alpha5_0 <- array(t(inits$alpha5_0), dim = c(tt, N/n_chains, n_chains))
  alpha5_0 <- aperm(alpha5_0, c(2,1,3))
  
  inits <- list(alpha0_inits = alpha0_0,
                alpha2_inits = alpha2_0,
                alpha1_inits = alpha1_0,
                alpha4_inits = alpha4_0,
                alpha5_inits = alpha5_0)
  
  RjF <- matrix(rep(rowSums(data_f_j), N/n_chains), nrow = N/n_chains)
  RaF <- matrix(rep(rowSums(data_f_a), N/n_chains), nrow = N/n_chains)
  RjM <- matrix(rep(rowSums(data_m_j), N/n_chains), nrow = N/n_chains)
  RaM <- matrix(rep(rowSums(data_m_a), N/n_chains), nrow = N/n_chains)
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA_f_j", "DATA_f_a", "DATA_m_j", "DATA_m_a", 
                      "inits", "model"),
                envir = environment())
  
  mcmc_list <- foreach(i=1:n_chains, .packages = c("rjags")) %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(MjF=DATA_f_j,
                                  MaF=DATA_f_a,
                                  MjM=DATA_m_j,
                                  MaM=DATA_m_a,
                                  RjF=rowSums(data_f_j),
                                  RaF=rowSums(data_f_a),
                                  RjM=rowSums(data_m_j),
                                  RaM=rowSums(data_m_a),
                                  tt=tt))
    
    jags <- coda.samples(model = foo,
                         variable.names = c("alpha0", "alpha2", "alpha1",
                                            "alpha4", "alpha5"),
                         n.iter = Ntotal)
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("alpha0", "alpha2", "alpha1", "alpha4", "alpha5")
  var_names_index <- 1:25
  
  var_names_full <- rep(NA, 29)
  var_names_full[1] <- paste0("^", var_names[1], "\\[") # alpha0
  var_names_full[2] <- paste0("^", var_names[2], "\\[") # alpha2
  var_names_full[3] <- paste0("^", var_names[3], "\\[") # alpha1
  var_names_full[4] <- paste0("^", var_names[4], "\\[") # alpha4
  for (i in 1:(tt-1)) {
    var_name_full[i+4] <- paste0("^", var_names[5], "\\[[0-9]+", var_names_index[i], "\\]") # alpha5
  }
  
  alpha0_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  alpha0_res <- unlist(alpha0_res)
  alpha0_res <- matrix(alpha0_res, nrow = Ntotal)[Ntotal,]
  
  alpha2_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[2], colnames(jags_out1))])
  alpha2_res <- unlist(alpha2_res)
  alpha2_res <- matrix(alpha2_res, nrow = Ntotal)[Ntotal,]
  
  alpha1_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[3], colnames(jags_out1))])
  alpha1_res <- unlist(alpha1_res)
  alpha1_res <- matrix(alpha1_res, nrow = Ntotal)[Ntotal,]
  
  alpha4_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[4], colnames(jags_out1))])
  alpha4_res <- unlist(alpha4_res)
  alpha4_res <- matrix(alpha4_res, nrow = Ntotal)[Ntotal,]
  
  alpha5_res <- matrix(rep(NA, (tt-1)*N), ncol = tt - 1)
  for (i in 1:(tt-1)) {
    alpha5_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+4], colnames(jags_out1))])
    alpha5_temp <- unlist(alpha5_temp)
    alpha5_res[,i] <- matrix(alpha5_temp, nrow = Ntotal)[Ntotal,]
  }
  
  return(list(alpha0_res = alpha0_res,
              alpha2_res = alpha2_res,
              alpha1_res = alpha1_res,
              alpha4_res = alpha4_res,
              alpha5_res = alpha5_res))
}











count_jags <- function(data,
                       alpha_j,
                       inits,
                       N,
                       Ntotal,
                       n_chains) {
  tt <- length(data)
  
  DATA <- matrix(rep(data, N/n_chains), nrow = N/n_chains, byrow = T)
  
  model <- "model{
    for (i in 1:50) {
      p[i] <- 1/50
    }
  
    for (i in 1:N) {
      logit(deltaJF[i]) <- alpha0[i]
      logit(deltaAF[i]) <- alpha0[i] + alpha2[i]  
  
      xJ[i,1] ~ dcat(p[1:50])
      sur[i,1] ~ dcat(p[1:50])
      imm[i,1] ~ dcat(p[1:50])
     
      alpha0[i] ~ dnorm(0, 1/4) T(-10, 10)
      alpha2[i] ~ dnorm(0, 1/4) T(-10, 10)
      alpha6[i] ~ dnorm(0, 1/4) T(-10, 10)
      eta_t[i] <- exp(alpha6[i])
       
      rho[i] ~ dunif(0, 10)
       
      for (tt in 2:t) {
        rateJ[i,tt-1] <- .5 * rho[i] * deltaJF[i] * x[i,tt-1]
        xJ[i,tt] ~ dpois(rateJ[i,tt-1])
        sur[i,tt] ~ dbin(deltaAF[i], x[i,tt-1])
        rate_imm[i,tt-1] <- x[i,tt-1] * eta_t[i]
        imm[i,tt] ~ dpois(rate_imm[i,tt-1])
      }
    
      for (tt in 1:t) {
        x[i,tt] <- xJ[i,tt] + sur[i,tt] + imm[i,tt]
        logL[i,tt] <- logdensity.pois(y[i,tt], x[i,tt])
        phi[i,tt] <- -alpha * logL[i,tt] + 1e20
        zeros[i,tt] ~ dpois(phi[i,tt])
      }
    }
  }"
  
  alpha0_0 <- matrix(inits$alpha0_0, ncol = n_chains)
  alpha2_0 <- matrix(inits$alpha2_0, ncol = n_chains)
  rho_0 <- matrix(inits$rho_0, ncol = n_chains)
  alpha6_0 <- matrix(inits$alpha6_0, ncol = n_chains)
  xJ_0 <- array(t(inits$xJ_0), dim = c(tt, N/n_chains, n_chains))
  xJ_0 <- aperm(xJ_0, c(2,1,3))
  sur_0 <- array(t(inits$sur_0), dim = c(tt, N/n_chains, n_chains))
  sur_0 <- aperm(sur_0, c(2,1,3))
  imm_0 <- array(t(inits$imm_0), dim = c(tt, N/n_chains, n_chains))
  imm_0 <- aperm(imm_0, c(2,1,3))
  
  inits <- list(alpha0_inits = alpha0_0,
                alpha2_inits = alpha2_0,
                rho_inits = rho_0,
                alpha6_inits = alpha6_0,
                xJ_inits = xJ_0,
                sur_inits = sur_0,
                imm_inits = imm_0)
  
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  clusterExport(cl, c("DATA", "inits", "model"),
                envir = environment())
  
  mcmc_list <- foreach(i=1:n_chains, .packages = c("rjags")) %dopar% {
    foo <- jags.model(textConnection(model),
                      data = list(y = DATA,
                                  t = tt,
                                  N = N/n_chains,
                                  alpha = alpha_j,
                                  zeros = matrix(rep(0, tt*N/n_chains), ncol=tt)),
                      inits = list(alpha0 = inits$alpha0_inits[,i],
                                   alpha2 = inits$alpha2_inits[,i],
                                   rho = inits$rho_inits[,i],
                                   alpha6 = inits$alpha6_inits[,i],
                                   xJ = inits$xJ_inits[,,i],
                                   sur = inits$sur_inits[,,i],
                                   imm = inits$imm_inits[,,i]))
    
    jags <- coda.samples(model = foo,
                         variable.names = c("alpha0", "alpha2", "rho",
                                            "alpha6", "xJ", "sur", "imm"),
                         n.iter = Ntotal)
  }
  
  stopCluster(cl)
  
  jags_out <- as.matrix(mcmc_list)
  jags_out1 <- as.matrix(mcmc_list[[1]])
  
  var_names <- c("alpha0", "alpha2", "rho", "alpha6", "xJ", "sur", "imm")
  var_names_index <- 1:26
  
  # Transform the variables names into the format that can be used by grep()
  # for partial matching
  var_names_full <- rep(NA, 82)
  var_names_full[1] <- paste0("^", var_names[1], "\\[") # alpha0
  var_names_full[2] <- paste0("^", var_names[2], "\\[") # alpha2
  var_names_full[3] <- paste0("^", var_names[3], "\\[") # rho
  var_names_full[4] <- paste0("^", var_names[4], "\\[") # alpha6
  for (i in 1:tt) {
    var_names_full[i+4] <- paste0("^", var_names[5], "\\[[0-9]+,", var_names_index[i], "\\]") # xJ
    var_names_full[i+tt+4] <- paste0("^", var_names[6], "\\[[0-9]+,", var_names_index[i], "\\]") # sur
    var_names_full[i+2*tt+4] <- paste0("^", var_names[7], "\\[[0-9]+,", var_names_index[i], "\\]") # imm
  }
  
  alpha0_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[1], colnames(jags_out1))])
  alpha0_res <- unlist(alpha0_res)
  alpha0_res <- matrix(alpha0_res, nrow = Ntotal)[Ntotal,]
  
  alpha2_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[2], colnames(jags_out1))])
  alpha2_res <- unlist(alpha2_res)
  alpha2_res <- matrix(alpha2_res, nrow = Ntotal)[Ntotal,]
  
  rho_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[3], colnames(jags_out1))])
  rho_res <- unlist(rho_res)
  rho_res <- matrix(rho_res, nrow = Ntotal)[Ntotal,]
  
  alpha6_res <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[4], colnames(jags_out1))])
  alpha6_res <- unlist(alpha6_res)
  alpha6_res <- matrix(alpha6_res, nrow = Ntotal)[Ntotal,]
  
  xJ_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    xJ_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+4], colnames(jags_out1))])
    xJ_temp <- unlist(xJ_temp)
    xJ_res[,i] <- matrix(xJ_temp, nrow = Ntotal)[Ntotal,]
  }
  
  sur_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    sur_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+tt+4], colnames(jags_out1))])
    sur_temp <- unlist(sur_temp)
    sur_res[,i] <- matrix(sur_temp, nrow = Ntotal)[Ntotal,]
  }
  
  imm_res <- matrix(rep(NA, tt*N), ncol = tt)
  for (i in 1:tt) {
    imm_temp <- lapply(mcmc_list, function(jags_out) jags_out[,grep(var_names_full[i+2*tt+4], colnames(jags_out1))])
    imm_temp <- unlist(imm_temp)
    imm_res [,i] <- matrix(imm_temp, nrow = Ntotal)[Ntotal,]
  }
  
  
  return(list(alpha0_res = alpha0_res,
              alpha2_res = alpha2_res,
              rho_res = rho_res,
              alpha6_res = alpha6_res,
              xJ_res = xJ_res,
              sur_res = sur_res,
              imm_res = imm_res))
}















