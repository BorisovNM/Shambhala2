library(stats)
library(rje)

   IFN = paste(prefix, "Q0.csv", sep = "")
   EXP = read.table(IFN, header = TRUE, sep = "," )
   
SYMBOL = as.vector(EXP[,1])
MAS0 = as.matrix(EXP[,-1])
NS = ncol(MAS0)
Q = EXP

IFN = paste(prefix, "SampleNames.csv", sep = "")
res = read.table(IFN, header = TRUE, sep = "," )
CN = as.vector(res[,1])

    IFN = "OutputP0.txt"
     
    MAS = read.table(IFN, header = FALSE, sep =" ")
    
    NS = ncol(MAS)

    MAS = as.matrix(MAS)

    MAS1 = merge(MAS,Q, by = 1)

    NS1 = ncol(MAS1)

    MAS2 = MAS1[,2:NS]

    NC2 = ncol(MAS2)
    NR2 = nrow(MAS2)

    MAS2 = as.matrix(MAS2)

    MAS20 = as.numeric(MAS2)

    MAS21 = matrix(MAS20, nrow = NR2, ncol = NC2)

    MAS3 =  MAS1[,(NS+1):NS1]

    RM = rowMeans(log(MAS3))
    RS = rowSds(log(MAS3)) 
    Rm = rowMins(log(MAS3))

    MAS22 = log(MAS21+1)

    NR = nrow(MAS22)
 
    for ( nr in 1:NR ) {

        MAS22[nr,] = RM[nr] + RS[nr]*MAS22[nr,]

    }

    MAS23 = exp(MAS22)

    SYMBOL = as.vector(MAS1[,1])

    OFN = paste(prefix, "Harmonized.csv", sep = "")
    write.table(MAS23, OFN, col.names = CN, row.names = SYMBOL, sep =",")

   
