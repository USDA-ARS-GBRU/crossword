haplotypes_out <- function(pop,output,parental_genotypes)
{
	if(missing(pop))
	{
		stop('ERROR: input haplotypes should be provided')
	}
	if(missing(output))
	{
		print('WARNING: output file prefex should be provided')
	}
	if(missing(parental_genotypes))
	{
		stop('ERROR: parental_genotypes should be provided')
	}
	founders = colnames(parental_genotypes[[1]])
	haplotypes = pop
	hap_ids = names(haplotypes$haplotype)
	chr_ids = names(haplotypes$haplotype[[1]])
	A=data.frame()
	count = 1
	for (x in 1:length(hap_ids))
	{
		for (y in 1:length(chr_ids))
		{
			a = get_haplotype_matrix(haplotypes$haplotype[[x]][[y]])	
			for (z in 1:dim(a)[1])
			{
				A[count,1] = hap_ids[x]
				A[count,2] = chr_ids[y]
				A[count,3] = a[z,1]
				A[count,4] = founders[ceiling(a[z,2]/2)]
				A[count,5] = founders[ceiling(a[z,3]/2)]
				count=count+1
			}	
		}
	}
	colnames(A) = c('haplotype_id','chromosome_id','genetic_loci','X1','X2')
	if (output !='')
	{
		write.table(A,paste0(output,'.csv'),col.names=TRUE,row.names=FALSE,quote=FALSE)
		dput(haplotypes,file=paste0(output,'.haplo'))
	}
}

