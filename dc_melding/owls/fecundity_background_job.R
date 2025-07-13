library(rjags)

fecundity <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/fecundity.dat")

source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/fecundity_model.R")

invisible(capture.output(out_fec <- fecundity_model(data = fecundity,
                                                     N = 8000,
                                                     Ntotal = 50)))

save(out_fec, file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/fecundity_result.RData")