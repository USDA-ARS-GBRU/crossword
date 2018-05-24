get_values <- function(geno,qtn_effect,dominant)
{
    if(missing(geno))
    {
        stop("ERROR: genotype object should be provided")
    }
    if(missing(qtn_effect))
    {
        stop("ERROR: qtn_effect object should be provided")
    }
    if (missing(dominant))
	{
	    dominant = FALSE
	}else
	{
	    dominant = dominant
	}    
    Value = NA
    for (y in 1:ncol(geno))
    {
        value = 0
        for (x in 1:nrow(qtn_effect))
        {
            if(as.character(geno[x,y]) == as.character(qtn_effect$ref[x]))
            {
                value=value+as.numeric(as.character(qtn_effect$effect[x]))
            }else if(mapply(adist,as.character(geno[x,y]),as.character(qtn_effect$ref[x])) == 2)
            {
                 value=value-as.numeric(as.character(qtn_effect$effect[x]))
            }else
            {
                if(dominant == TRUE)
                {
                    value=value+as.numeric(as.character(qtn_effect$effect[x]))
                }else
                {
                    value=value+0
                }
            }
        }
        Value=c(Value,value)
    }
    Value <- Value[!is.na(Value)]
    names(Value) = colnames(geno)
    return(Value)
}

