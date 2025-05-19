
cap_recap_model <- function(data_f_j,
                            data_f_a,
                            data_m_j,
                            data_m_a,
                            N,
                            Ntotal) {
  
  time_begin <- Sys.time()
  
  model <- "model {
    
    #----------Part 1. Regression equations (44)-------------#
    ##------------------------------------------------------------##
    ## alpha = (alpha0, alpha1, alpha2, alpha4)
    ##------------------------------------------------------------##
    for (i in 1:(ti-1)) {
      # Recapture rate pi
      logit(pM[i]) <- alpha[4] + alpha5[i]   # Male
      logit(pF[i]) <- alpha5[i]              # Female
    }
    
    # Juvenile survival rate deltaj
    logit(deltajM) <- alpha[1] + alpha[2]  # Male
    logit(deltajF) <- alpha[1]             # Female  
    
    # Adult survival rate deltaa
    logit(deltaaM) <- alpha[1] + alpha[2] + alpha[3]   # Male
    logit(deltaaF) <- alpha[1] + alpha[3]              # Female
    
    
    #----------Part 2. Priors on alpha's-------------#
    for (i in 1:4) {
      alpha[i] ~ dnorm(0, 1/4) T(-10, 10)
    }
    
    for (i in 1:(ti-1)) {
      alpha5[i] ~ dnorm(1, 1/4) T(-10, 10)
    }
    
    
    #----------Part 3. Likelihood-------------#
    ##------------------------------------------------------------##
    ## Female capture recapture rate
    ##------------------------------------------------------------##
    ### Rate for juveniles ###
    for (i in 1:(ti-1)) {
      #RjF[i] <- sum(MjF[i, 1:ti])
      MjF[i, 1:ti] ~ dmulti(QjF[i, 1:ti], RjF[i])
      qF[i] <- 1 - pF[i]  # 1 - pi
    }
    
    #### Definition on QjF ####
    for (i in 1:(ti-1)) {
      ###qF[i] <- 1 - pF[i]  # 1 - pi -- or it should be qF <- 1 - pF outside the loop???
      
      # Main diagonal
      QjF[i, i] <- deltajF * pF[i]
      
      # Above main diagonal
      for (j in (i+1):(ti-1)) {
        QjF[i, j] <- deltajF * deltaaF * pF[j] * prod(qF[i:(j-1)]) # Is deltajF^j or just deltajF?
      }
      
      # Below main diagonal
      for (j in 1:(i-1)) {
        QjF[i, j] <- 0
      }
      
      # Last column
      QjF[i, ti] <- 1 - sum(QjF[i, 1:(ti-1)])
    }
    
    
    ### Rate for adults ###
    for (i in 1:(ti-1)) {
      #RaF[i] <- sum(MaF[i, 1:ti])
      MaF[i, 1:ti] ~ dmulti(QaF[i, 1:ti], RaF[i])
    }
    
    #### Definition on QaF ####
    for (i in 1:(ti-1)) {
      ###qF[i] <- 1 - pF[i]  # Do I need define here again?
      
      # Main diagonal
      QaF[i, i] <- deltaaF * pF[i]
      
      # Above main diagonal
      for (j in (i+1):(ti-1)) {
        QaF[i, j] <- deltaaF * pF[j] * prod(qF[i:(j-1)]) # Is deltaaF^j or just deltaaF?
      }
      
      # Below main diagonal
      for (j in 1:(i-1)) {
        QaF[i, j] <- 0
      }
      
      # Last column
      QaF[i, ti] <- 1 - sum(QaF[i, 1:(ti-1)])
    }
    
    
    
    ##------------------------------------------------------------##
    ## Male capture recapture rate
    ##------------------------------------------------------------##
    ### Rate for juveniles ###
    for (i in 1:(ti-1)) {
      #RjM[i] <- sum(MjM[i, 1:ti])
      MjM[i, 1:ti] ~ dmulti(QjM[i, 1:ti], RjM[i])
      qM[i] <- 1 - pM[i]
    }
    
    #### Definition on QjM ####
    for (i in 1:(ti-1)) {
      ###qM[i] <- 1 - pM[i]  # 1 - pi -- or it should be qM = 1 - pM outside the loop???
      
      # Main diagonal
      QjM[i, i] <- deltajM * pM[i]
      
      # Above main diagonal
      for (j in (i+1):(ti-1)) {
        QjM[i, j] <- deltajM * deltaaM * pM[j] * prod(qM[i:(j-1)]) # Is deltajM^j or just deltajF?
      }
      
      # Below main diagonal
      for (j in 1:(i-1)) {
        QjM[i, j] <- 0
      }
      
      # Last column
      QjM[i, ti] <- 1 - sum(QjM[i, 1:(ti-1)])
    }
    
    
    ### Rate for adults ###
    for (i in 1:(ti-1)) {
      #RaM[i] <- sum(MaM[i, 1:ti])
      MaM[i, 1:ti] ~ dmulti(QaM[i, 1:ti], RaM[i])
    }
    
    #### Definition on QaF ####
    for (i in 1:(ti-1)) {
      ###qM[i] <- 1 - pM[i]  # Do I need define here again?
      
      # Main diagonal
      QaM[i, i] <- deltaaM * pM[i]
      
      # Above main diagonal
      for (j in (i+1):(ti-1)) {
        QaM[i, j] <- deltaaM * pM[j] * prod(qM[i:(j-1)]) # Is deltaaM^j or just deltaaF?
      }
      
      # Below main diagonal
      for (j in 1:(i-1)) {
        QaM[i, j] <- 0
      }
      
      # Last column
      QaM[i, ti] <- 1 - sum(QaM[i, 1:(ti-1)])
    }
    
  }"
  
  
  parameters <- c("alpha")
  
  foo <- jags.model(textConnection(model),
                    data = list(MjF=data_f_j,
                                MaF=data_f_a,
                                MjM=data_m_j,
                                MaM=data_m_a,
                                RjF=rowSums(data_f_j),
                                RaF=rowSums(data_f_a),
                                RjM=rowSums(data_m_j),
                                RaM=rowSums(data_m_a),
                                ti=dim(data_f_j)[2]),
                    quiet = T)
  update(foo)
  
  alpha <- array(rep(NA, Ntotal*4*N), dim = c(Ntotal, 4, N))
  #alpha5 <- array(rep(NA, N*Ntotal), dim = c(N, Ntotal))
  
  for (i in 1:N) {
    jags_res <- coda.samples(model = foo, variable.names = parameters,
                             n.iter = Ntotal, thin = 1)
    out_jags <- as.matrix(jags_res)
    alpha[,,i] <- out_jags
  }
  
  time_end <- Sys.time()
  time_run <- time_end - time_begin
  
  return(list(alpha=alpha, tim=time_run))
}



