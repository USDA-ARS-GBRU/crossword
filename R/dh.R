dh <- function(pop)
{
    A = pop[[1]]
    gen = A$gen[nrow(A)]
    B = pop[[2]]
    CHR = names(pop[[2]][[1]])
    for (i in 1:nrow(A))
    {
        if(A$gen[i]==gen)
        {
            id1 = A[i,1]
            for (y in 1:length(CHR))
            {
                B[[id1]][[CHR[y]]][[2]] = B[[id1]][[CHR[y]]][[1]]
            }
        }
    }
    return(list(scheme=A,haplotype=B))
}

