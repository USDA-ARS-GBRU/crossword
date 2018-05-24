combine_populations <- function(pops)
{
	if(missing(pops))
	{
		stop('population objects should be provided')
	}
	if(length(pops) <= 1)
	{
	    stop("ERROR, populations should be more than one")
	}
	scheme = pops[[1]]
	haplotype = pops[[2]]
	for(i in seq(1,length(pops),by=2))
	{
	    if(i == 1){next}
	    scheme = rbind(scheme,pops[[i]])
	    haplotype = c(haplotype,pops[[i+1]])
	}
	Pops = list(scheme=scheme,haplotype=haplotype)
	return(Pops)
}

