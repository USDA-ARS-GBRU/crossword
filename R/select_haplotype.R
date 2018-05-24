select_haplotype <- function(pop,haplotypes_ids)
{
	if(missing(pop))
	{
		stop('ERROR: haplotypes should be provided')
	}
	if(missing(haplotypes_ids))
	{
		stop('ERROR: The IDs of haplotypes should be provided')
	}else if(is.numeric(haplotypes_ids))
	{
	    hids = pop[[1]][haplotypes_ids,1]
	}else
	{
	    hids = haplotypes_ids
	}
	haplotype = pop
	A = data.frame()
	B = list()
	for (i in hids)
	{
		a = haplotype[[1]][haplotype[[1]]$id==i,]
		b = haplotype[[2]][[i]]
		A[i,1] = a[1]
		A[i,2] = a[2]
		A[i,3] = a[3]
		A[i,4] = a[4]
		A[i,5] = a[5]
		A[i,6] = a[6]
		A[i,7] = a[7]
		B[[i]] = b
	}
	return(list(scheme = as.data.frame(A), haplotype = B))
}

