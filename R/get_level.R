get_level <- function(m1,level,return_list,print_warnings)
{  
    ### check if the return will be one dimension array or a list
    if(missing(return_list))
    {
        return_list = FALSE
    }else
    {
        return_list = return_list
    }
    if(missing(print_warnings))
    {
        print_warnings = FALSE 
    }else
    {
        print_warnings = print_warnings
    }
    if(level=="individual")
	{
	    l1 = paste0(as.character(m1$gen),"_",as.character(m1$C),"_",as.character(m1$S),"_",as.character(m1$N))
	}else if(level=="family")
	{
	    l1 = paste0(as.character(m1$gen),"_",as.character(m1$C),"_",as.character(m1$S))
	    if (sum(is.na(m1$S) ) > 0)
	    {
	        stop("ERROR, you selected level = family and you do not have a family structure")
	    }
	}else if(level=="cross")
	{
	    l1 = paste0(as.character(m1$gen),"_",as.character(m1$C))
	    if(length(unique(l1)) == 1)
	    {
	        if(print_warnings == TRUE)
	        {
	            print("WARNING: you selected level = cross, and you only has one cross in this population")
	        }
	    }
	}else if(level=="population")
	{
	    l1 = as.character(m1$gen)
	}else 
	{
	    stop("ERROR, you selected an incorrect level")
	}
	l1 = gsub(" ","",l1)
	if(return_list == TRUE)
	{
	    L1 = list()
	    for (x in unique(l1))
        {
            L1[[x]] = m1[l1 == x,]
        }
        return(L1)
    }else
    {
        return(l1)
    }
}

