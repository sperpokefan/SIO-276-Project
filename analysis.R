ETH <- read.delim(file = "ETH.dat", header = FALSE)
LTC <- read.delim(file = "LTC.dat", header = FALSE)
XMR <- read.delim(file = "XMR.dat", header = FALSE)
XRP <- read.delim(file = "XRP.dat", header = FALSE)

ETH_vect <- ETH$V1[seq(1, length(ETH$V1)/30, 6)]
LTC_vect <- LTC$V1[seq(1, length(LTC$V1)/30, 6)]
XMR_vect <- XMR$V1[seq(1, length(XMR$V1)/30, 6)]
XRP_vect <- XRP$V1[seq(1, length(XRP$V1)/30, 6)]

ts_matrix <- as.data.frame(matrix(NA, 96, 5))

colnames(ts_matrix) <- c("Time", "ETH", "LTC", "XMR", "XRP")

ts_matrix$Time <- c(1:96)
ts_matrix$ETH <- ETH_vect
ts_matrix$LTC <- LTC_vect
ts_matrix$XMR <- XMR_vect
ts_matrix$XRP <- XRP_vect


lib <- c(1,96)
pred <- lib

save(ts_matrix, lib, pred, file = "crypto_ts_matrix.Rdata")

source("simplex_smap.R")


simplex_smap_bestE(in_file = "crypto_ts_matrix.Rdata", out_file = "crypto_simplex_smap.Rdata", Best_E_File = "crypto_best_E.csv")

source("ccm_network.R")

compute_ccm_web(in_file = "crypto_ts_matrix.Rdata", best_E_file = "crypto_best_E.csv", num_surr = 1000, out_file = "crypto_ccm_network.Rdata")


load("crypto_ccm_network.Rdata")

