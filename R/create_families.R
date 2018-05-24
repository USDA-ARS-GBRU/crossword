create_families <- function(pop,S,id,chr_length,is_clone)
{
    if(missing(is_clone) || is.na(is_clone))
    {
        is_clone = FALSE
    }else
    {
        is_clone = is_clone
    }
    C0 = pop[[1]]
    C1 = C0[C0$gen == C0[nrow(C0),2],]
    A = matrix(nrow=S*nrow(C1),ncol=7)
    count_x = 1
    gen = as.numeric(as.character(C0[nrow(C0),2])) + 1
    for (i in 1:nrow(C1))
    {
        for (i2 in 1:S)
        {
            A[count_x,1] = gsub(" ", "",paste0(id,"_",count_x,"_",as.character(C1[i,5]),"_",as.character(C1[i,7]),"_",i2))
            A[count_x,2] = gen
            A[count_x,3:6] = c(as.character(C1[i,1]),as.character(C1[i,1]),C1[i,5],C1[i,7])
            A[count_x,7] = i2
            count_x = count_x+1
        }
    }
    colnames(A) = colnames(C0)
    A = as.data.frame(A)
    B = list()
    if(is_clone == FALSE)
    {
        haplotype = create_haplotypes_from_scheme(A,chr_length,pop)
	    return(list(scheme=A,haplotype=haplotype))
	}else
	{
	    for (i in 1:nrow(A))
	    {
	        haplo = select_haplotype(pop,as.character(A[i,3]))
	        haplo[[1]] = A[i,]
	        names(haplo[[2]])[1] = as.character(A$id[i])
	        if(i == 1)
	        {
	            B = haplo
	        }else
	        {
	            B = combine_populations(c(B,haplo))
	        }
	    }
	    return(B)
	}
}

