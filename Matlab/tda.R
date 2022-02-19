rm(list = ls())

library(R.matlab)

data <- readMat("data/MPRC_transformed_cor.mat")
print(data)

X_N = data$cor.NC.fisherZ.correct
X_S = data$cor.SZ.fisherZ.correct

library(TDA) 

X_N_pd <- list()
for(i in seq(1,dim(X_N)[3])){
  print(i)
  X_N_pd[i] = gridDiag(FUNvalues = X_N[,,i])
  # plot(X_N_pd[[i]])
}

X_S_pd <- list()
for(i in seq(1,dim(X_S)[3])){
  print(i)
  X_S_pd[i] = gridDiag(FUNvalues = X_S[,,i])
  # plot(X_S_pd[[i]])
}

writeMat("data/MPRC_transformed_cor_persistence_diagram.mat", X_N_pd=X_N_pd, X_S_pd=X_S_pd)

