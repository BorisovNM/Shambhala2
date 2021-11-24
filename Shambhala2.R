library(matrixStats)

# k (OPTIONAL) is the number of probe clusters for the application of k-means to find probe-cluster partitions. By default it is 5.

Shambhala2 <- function( InputFileName, OutputFileName, PFileName, QFileName, delete_buffer_files = TRUE, k = 5 ) {

    IFN = InputFileName
    MAS = read.table(IFN, header = TRUE, sep = ",")
    MAS0 = as.matrix(MAS[,-1])
    NS = ncol(MAS0)
    SYMBOL = as.vector(MAS[,1])
    CN = colnames(MAS)

    PFN = PFileName 
    P = read.table(PFN, header = TRUE, sep = ",")
   
    pool = merge(MAS,P,by = "SYMBOL")

    NS1 = ncol(pool)

    NG1 = nrow(pool)

    for ( j in (NS+2):NS1 ) {
        pool[,j] = as.numeric(pool[,j])
    }     
   
    P1FN = "P_prim.txt"
    write.table(pool, P1FN, row.names = FALSE, col.names = TRUE, sep = "\t")

    NH = ncol(MAS) - 1
    NP = ncol(P) - 1

    args = c(NH,NP,k)
    
    AFN = "args.txt"
    write.table(args, AFN, row.names = FALSE, col.names = FALSE)

    system("matlab -nodesktop -nosplash -nodisplay -r \"run('Shambhala2.m');exit;\"")

    QFN = QFileName
    Q = read.table(QFN, header = TRUE, sep = ",")

    Cu2FN = "Cu_bis.txt"
    MAS = read.table(Cu2FN, header = FALSE, sep = " ")

    MAS = as.matrix(MAS)

    MAS1 = merge(MAS,Q, by = 1)

    NS1 = ncol(MAS1)

    for ( j in 2:NS1 ) {
        MAS1[,j] = as.numeric(MAS1[,j])
    } 

    MAS3 = MAS1[,(NS+2):ncol(MAS1)]

    RM = rowMeans(log(as.matrix(MAS3)))
    RS = rowSds(log(as.matrix(MAS3))) 

    MAS2 = MAS1[,2:(NS+1)]

    MAS22 = log(MAS2+1)

    NR = nrow(MAS22)
 
    for ( nr in 1:NR ) {
        MAS22[nr,] = RM[nr] + RS[nr]*MAS22[nr,]
    }

    MAS23 = exp(MAS22)

    SYMBOL = as.vector(MAS1[,1])

    MAS33 = cbind(SYMBOL,MAS23)

    for ( j in 2:(NS+1) ) {
        MAS33[,j] = as.vector(as.numeric(MAS33[,j]))
    } 

    OFN = OutputFileName
    write.table(as.matrix(MAS33), OFN, col.names = CN, row.names = FALSE, sep =",")

    if ( delete_buffer_files ) {

        if (file.exists(P1FN)) {
        #Delete file if it exists
            file.remove(P1FN) 
        }
    
        if (file.exists(Cu2FN)) {
        #Delete file if it exists
            file.remove(Cu2FN) 
        }
  
        if (file.exists(AFN)) {
        #Delete file if it exists
            file.remove(AFN) 
        }

    }
    
    colnames(MAS33) = CN
    
    return(MAS33)
    
}    
   
Harmonized = Shambhala2("Input.csv", "Output.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5) 
