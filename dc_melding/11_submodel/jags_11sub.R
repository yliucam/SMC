library(rjags)


model <- "model{
  phi12 ~ dnorm(0, 1/4)
  phi23 ~ dunif(1e-5, 30) #dgamma(2, 2)
  phi34 ~ dunif(-100, 100) #dnorm(0, 1/4)
  phi45 ~ dgamma(2, 2)
  
  psi1 ~ dunif(1e-5, 50) #dgamma(2, 2)
  psi2 ~ dcat(nu_p[])
  psi3 ~ dcat(nu_p[])
  psi4 ~ dbeta(9, 1)
  
  phi10_11 ~ dgamma(2, 2)
  phi9_10 ~ dunif(-100, 100) #dnorm(0, 1/4)
  phi89 ~ dunif(1e-5, 30) #dgamma(2, 2)
  phi78 ~ dnorm(0, 1/4)
  
  psi11 ~ dnorm(0, 1/4)
  psi10 ~ dcat(nu_p[])
  psi9 ~ dcat(nu_p[])
  psi8 ~ dbeta(9, 1)
  
  phi56 ~ dgamma(2, 2)
  phi67 ~ dnorm(0, 1/4)
  
  psi6 ~ dbeta(9, 1)
  
  # model 1, 2 and 3
  for (i in 1:n) {
    y1[i] ~ dnorm(phi12, 1/psi1^2)
    y2[i] ~ dt(phi12, phi23, psi2)
    y3[i] ~ dt(phi34, phi23, psi3)
  }
  
  # model 4
  tau45 <- pow(phi45, -2)
  x4_0 ~ dnorm(phi34, tau45)
  x4mean[1] <- phi34 + psi4 * (x4_0 - phi34)
  x4[1] ~ dnorm(x4mean[1], tau45)
  y4[1] ~ dnorm(x4[1], 1)
  for (t in 2:TT) {
    x4mean[t] <- phi34 + psi4 * (x4[t-1] - phi34)
    x4[t] ~ dnorm(x4mean[t], tau45)
    y4[t] ~ dnorm(x4[t], 1)
  }
  
  # model 5
  tau56 <- pow(phi56, -2)
  x5_0 ~ dnorm(1, tau56)
  x5[1] ~ dnorm(x5_0, tau56)
  y5[1] ~ dnorm(x5[1], tau45)
  for (t in 2:TT) {
    x5[t] ~ dnorm(x5[t-1], tau56)
    y5[t] ~ dnorm(x5[t], tau45)
  }
  
  # model 9, 10 and 11
  for (i in 1:n) {
    y11[i] ~ dnorm(psi11, 1/phi10_11^2)
    y10[i] ~ dt(phi9_10, phi10_11, psi10)
    y9[i] ~ dt(phi9_10, phi89, psi9)
  }
  
  # model 8
  tau89 <- pow(phi89, -2)
  x8_0 ~ dnorm(phi78, tau89)
  x8mean[1] <- phi78 + psi8 * (x8_0 - phi78)
  x8[1] ~ dnorm(x8mean[1], tau89)
  y8[1] ~ dnorm(x8[1], 1)
  for (t in 2:TT) {
    x8mean[t] <- phi78 + psi8 * (x8[t-1] - phi78)
    x8[t] ~ dnorm(x8mean[t], tau89)
    y8[t] ~ dnorm(x8[t], 1)
  }
  
  # model 7
  x7_0 ~ dlnorm(phi67, 1/.01)
  x7[1] ~ dlnorm(x7_0 + phi67, 1/.01)
  tau7[1] <- pow(x7[1], -2) + 1e-10
  y7[1] ~ dnorm(phi78, tau7[1])
  for (t in 2:TT) {
    x7[t] ~ dlnorm(x7[t-1] + phi67, 1/.01)
    tau7[t] <- pow(x7[t], -2) + 1e-10
    y7[t] ~ dnorm(phi78, tau7[t])
  }
  
  # model 6
  x6_0 ~ dnorm(phi67, tau56)
  x6mean[1] <- phi67 + psi6 * (x6_0 - phi67)
  x6[1] ~ dnorm(x6mean[1], tau56)
  y6tau[1] <- 1/exp(x6[1])
  y6[1] ~ dnorm(0, y6tau[1])
  for (t in 2:TT) {
    x6mean[t] <- phi67 + psi6 * (x6[t-1] - phi67)
    x6[t] ~ dnorm(x6mean[t], tau56)
    y6tau[t] <- 1/exp(x6[t])
    y6[t] ~ dnorm(0, y6tau[t])
  }
}"


data_list <- list(y1 = y1, y2 = y2, y3 = y3, y4 = y4, 
                  y5 = y5, y6 = y6, y7 = y7, y8 = y8, 
                  y9 = y9, y10 = y10, y11 = y11,
                  n = 50, TT = 10, nu_p = rep(1/30, 30))

parameter <- c("phi12", "phi23", "phi34", "phi45", "phi56", "phi67",
               "phi78", "phi89", "phi9_10", "phi10_11",
               "psi1", "psi2", "psi3", "psi4", "psi6",
               "psi8", "psi9", "psi10", "psi11")

foo <- jags.model(textConnection(model), 
                  data = data_list,
                  inits = list(x7 = rep(1e-5, 10)))

update(foo, 5000)

out_jags <- coda.samples(foo,
                         variable.names = parameter,
                         n.iter = 10000)
out_jags <- as.matrix(out_jags)





