library(rjags)

recap_f_a <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/capRecapFemaleAdult.dat")
recap_f_j <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/capRecapFemaleFirst.dat")
recap_m_a <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/capRecapMaleAdult.dat")
recap_m_j <- read.table("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/data/capRecapMaleFirst.dat")

source("C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/cap_recap_model.R")

invisible(capture.output(out_recap <- cap_recap_model(data_f_j = recap_f_j, data_f_a = recap_f_a, 
                                                data_m_j = recap_m_j, data_m_a = recap_m_a, 
                                                N=8000, Ntotal = 50)
))

save(out_recap, file = "C:/Users/Yixuan/Documents/codes/SMC/dc_melding/owls/recap_result.RData")