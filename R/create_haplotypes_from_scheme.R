create_haplotypes_from_scheme <- function(scheme,chr_length,pop)
{
    require('simcross')
    if(missing(pop))
    {
    	haplotype=list()
    }else
    {
    	haplotype=pop$haplotype
    }
    CHR = read.table(chr_length)
    colnames(CHR) = c('id','legnth')
	rownames(CHR) = CHR$id
	for (x in 1:nrow(scheme))
	{
		for(y in 1:nrow(CHR))
		{
			id1 = as.character(scheme[x,1])
			if(id1 %in% pop[[1]]$id)
			{
			    haplotype[[id1]][[rownames(CHR[y,])]] = pop$haplotype[[id1]][[rownames(CHR[y,])]]

			}else
			{
			    haplotype[[id1]][[rownames(CHR[y,])]] = simcross::cross((haplotype[[as.character(scheme[x,3])]][[rownames(CHR[y,])]]), (haplotype[[as.character(scheme[x,4])]][[rownames(CHR[y,])]]))
			}
		}
	}
	return(haplotype[as.character(scheme$id)])
}

