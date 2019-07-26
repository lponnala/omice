library(readxl)

D = read_excel("SPC-PREP-OOP-forLalitv1.xlsx", sheet=1, skip=1)
D = D[-tail(seq(1,nrow(D)),2),]

Acc = D[,2]
save(Acc, file="Acc.RData")

N = D[,seq(30,41)]
names(N) = c("rep1-wt","rep1-prep","rep1-oop","rep1-triple",
						"rep2-wt","rep2-prep","rep2-oop","rep2-triple",
						"rep3-wt","rep3-prep","rep3-oop","rep3-triple")
save(N, file="N.RData")

totalA = D[,17]
save(totalA, file="totalA.RData")

A = D[,seq(18,29)]
names(A) = c("rep1-wt","rep1-prep","rep1-oop","rep1-triple",
						"rep2-wt","rep2-prep","rep2-oop","rep2-triple",
						"rep3-wt","rep3-prep","rep3-oop","rep3-triple")
save(A, file="A.RData")


