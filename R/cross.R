## the function is used to create haplotype information based on the input
#the output is a list structure of the original scheme and the haplotype information
cross <- function (P1,P2,N,id,chr_length,parental_genotypes,isdh)
{
	require('simcross')
	if(is.list(P1))
	{
        PP1 = P1$scheme[P1$scheme$gen==P1$scheme$gen[nrow(P1$scheme)],]
        ids1 = as.character(PP1[,1])
		POP1 = select_haplotype(P1,ids1)
		p1 = ids1
		
	}else{p1=P1}
	if(missing(P2) || missing(P1))
	{
		stop("ERROR: the two parents should be provided")
	}
	if(is.list(P2))
	{
		PP2 = P2$scheme[P2$scheme$gen==P2$scheme$gen[nrow(P2$scheme)],]
		ids2 = as.character(PP2[,1])
		POP2 = select_haplotype(P2,ids2)
		p2 = ids2

	}else{p2=P2}
	##p1 and p2 holds the IDs of the recurrnt parents
	##P1 and P2 holds the haplotypes of the recurrent parents
	scheme = create_scheme(p1,p2,N,id)
	if(missing(isdh))
	{
	    isDH = FALSE
	}else
	{
	    isDH = isdh
	}
    founders = colnames(parental_genotypes[[1]])
    CHR = read.table(chr_length)
	colnames(CHR) = c('id','legnth')
	rownames(CHR) = CHR$id
	founders = data.frame(names=founders, id=seq(1,2*(length(founders)),by=2))
	rownames(founders) = founders[,1]	
	haplotype=list()
	for (x in 1:nrow(scheme))
	{
		for(y in 1:nrow(CHR))
		{
			id1 = as.character(scheme[x,1])
			if(is.list(P1) && id1 %in% names(POP1$haplotype))
			{
				haplotype[[id1]][[rownames(CHR[y,])]] = POP1$haplotype[[which(id1 %in% names(POP1$haplotype))]][[rownames(CHR[y,])]]
			}else if(is.list(P2) && id1 %in% names(POP2$haplotype))
			{
				haplotype[[id1]][[rownames(CHR[y,])]] = POP2$haplotype[[which(id1 %in% names(POP2$haplotype))]][[rownames(CHR[y,])]]
			}else
			{
				if(scheme[x,3][[1]] == 0 && scheme[x,4][[1]] == 0)
				{
					haplotype[[id1]][[rownames(CHR[y,])]] = create_parent(CHR[y,2], (founders[which(as.character(scheme[x,1]) == rownames(founders)),2]):((founders[which(as.character(scheme[x,1]) == rownames(founders)),2])+1))	
				}else
				{
					haplotype[[id1]][[rownames(CHR[y,])]] = simcross::cross((haplotype[[as.character(scheme[x,3])]][[rownames(CHR[y,])]]), (haplotype[[as.character(scheme[x,4])]][[rownames(CHR[y,])]]),obligate_chiasma=TRUE)
					if (isDH==TRUE)
					{
					    haplotype[[id1]][[rownames(CHR[y,])]][[1]] = haplotype[[id1]][[rownames(CHR[y,])]][[2]]
					}
				}
			}
		}
	}
	return(list(scheme=scheme,haplotype=haplotype))
}

