create_population <- function(P,chr_length,parental_genotypes)
{
	if(missing(P))
	{
	    stop("ERROR: P object should be provided")
	}
	if(missing(chr_length))
	{
	    stop("ERROR: chr_length file should be provided")
	}
	if(missing(parental_genotypes))
	{
	    stop("ERROR: parental_genotypes object should be provided")
	}
	require('simcross')
    founders = colnames(parental_genotypes[[1]])
    CHR = read.table(chr_length)
	colnames(CHR) = c('id','legnth')
	rownames(CHR) = CHR$id
	founders = data.frame(names=founders, id=seq(1,2*(length(founders)),by=2))
	rownames(founders) = founders[,1]	
	haplotype=list()
	scheme = data.frame(matrix(nrow=length(P), ncol=7))
	colnames(scheme) = c("id","gen","p1","p2","C","N","S")
	scheme[,1] = P
	scheme[,2:4] = 0
	for (x in 1:length(P))
	{
		for(y in 1:(dim(CHR)[1]))
		{
			    id1 = scheme[x,1][[1]]
				if(scheme[x,3][[1]] == 0 && scheme[x,4][[1]] == 0)
				{
					haplotype[[id1]][[rownames(CHR[y,])]] = create_parent(CHR[y,2], (founders[which(scheme[x,1][[1]] == rownames(founders)),2]):((founders[which(scheme[x,1][[1]] == rownames(founders)),2])+1))
					
				}
		}
	}
	return(list(scheme=scheme,haplotype=haplotype))
}

