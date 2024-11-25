library(splatter)
rate=-0.5
mid=3 #3 , 4 or 5
params = newSplatParams ( )
params = setParams(params,list(batchCells = 1000,
                               nGenes = 2000,
                               group.prob = rep(0.2 ,5),
                               de.prob = c (0.05,0.08,0.01,0.1,0.1),
                               de.facLoc = 0.5,
                               de.facScale = 0.8,
                               seed= 1))
sim = splatSimulateGroups(params,
                          dropout.shape = rep(rate,5) ,
                          dropout.mid = rep(mid,5) ,
                          dropout.type = "group" )

counts=as.data.frame(as.matrix(counts(sim)))
truecounts=as.data.frame(as.matrix(assays(sim)$TrueCounts))
dropout=as.data.frame(as.matrix(assays(sim)$Dropout))

write.csv(counts,"dataset/sim1/counts.csv")
write.csv(truecounts,"dataset/sim1/truecounts.csv")
write.csv(dropout,"dataset/sim1/dropout.csv")