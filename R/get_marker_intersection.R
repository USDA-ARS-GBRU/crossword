get_marker_intersection <- function (genotypes,input_geno_file,pheno,geno)
{
    if(missing(genotypes))
    {
        stop("ERROR: input genotypes should be provided")
    }
    A = genotypes[[1]]
    B = genotypes[[2]]
    if(!missing(input_geno_file))
    {
        A3 = read.table(input_geno_file,header=FALSE)
        C = paste0(A3[,1],"_",A3[,2])
    }else if (!missing(pheno))
    {
        C = rownames(pheno$affected_genotypes)
    }else if (!missing(geno))
    {
        C = rownames(geno[[1]])
    }
    A2 = A[C,]
    B2 = as.data.frame(B[C,])
    rownames(B2) = C
    colnames(B2) = "gen_loc"
    return(list(genotypes=A2,gen2phy=B2))
}

