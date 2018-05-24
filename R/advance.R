advance <- function(pop,F,outcross,level,clevel,id,chr_length)
{
    if(missing(pop))
	{
		stop('ERROR: a cross population s
		hould be given')
	}
	if(missing(F))
	{
		stop('ERROR: Number of generation should be provided')
	}	
	if(missing(outcross))
	{
		selfing = 1
		outcross= 0
	}else
	{
		selfing = 1-outcross
	}
	if(missing(level))
	{
		stop('ERROR: the level should be given')
	}
	if(missing(clevel))
	{
		stop('ERROR: the crossing strtegy should be given')
	}	
	####getting the scheme to be advanced
	C0 = pop$scheme
    count_x = 1
    for (y in 2:F)
    {
    	C1 = C0[C0$gen==C0[nrow(C0),2],]
    	L1 = get_level(C1,level,TRUE)
    	CL1 = get_level(C1,clevel,TRUE) 
	    if(y != F)
	    {
	        l1 = get_level(C1,level,FALSE)
	        cl1 = get_level(C1,clevel,FALSE)
	    }else
	    {
	        l1 = get_level(C1,level,FALSE,TRUE)
	        cl1 = get_level(C1,clevel,FALSE,TRUE)
	    }
	    C2 = as.matrix(C1) 
	    for (x in 1:nrow(C1))
	    {
	        mother_ind = as.character(sample(L1[[l1[x]]]$id,1))
	        if (sample(1:2,1,prob=c(selfing,outcross))==1)
		    {
		        father_ind = mother_ind
		    }else
		    {
		        father_ind = as.character(sample(CL1[[cl1[x]]]$id,1))
		    }
		    ID = paste0(id,"_",count_x,"_",sub(".*_.*_(.*_.*_.*)","\\1",mother_ind))
		    C2[x,1] = ID
		    C2[x,2] = y
		    C2[x,3] = mother_ind
		    C2[x,4] = father_ind
		    count_x = count_x +1
	    }
	    C0 = rbind(C0,as.data.frame(C2))
    }
    C0 = C0[(nrow(pop$scheme)+1):nrow(C0),]
    ####reading scheme to get individuals
    haplotype = create_haplotypes_from_scheme(C0,chr_length,pop)
    scheme = as.data.frame(C2)
    haplotype = haplotype[as.character(scheme$id)]
	return(list(scheme=scheme,haplotype=haplotype))
}

