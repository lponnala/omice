
load("N.RData")

write.csv(cor(N), file="pairwise-corr.csv")
