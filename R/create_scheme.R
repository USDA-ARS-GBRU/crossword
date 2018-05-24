### The function creates scheme from the provided information:
## 1. single parent1 X single parent2
## 2. single parent1 X list of parents2
## 3. list of parents1 X single parent2
# N: number of F1 crosses
create_scheme <- function(P1,P2,N,id)
{
	if(missing(P1))
	{
		stop('ERROR: parent/parents 1 should be given')
	}
	if(missing(P2))
	{
		stop('ERROR: parent/parents 2 should be given')	
	}
	if(missing(N))
	{
		stop('ERROR: Number of crosses per generation should be provided')
	}
	if(length(P1) > 1 && length(P2) > 1 )
	{
	    stop("ERROR: you can not cross population by poulation")
	}
    ###########geting parents
	parents = c(P1,P2)
	nparents = length(parents)
	A = matrix(nrow=nparents,ncol=7)
	A[1:nparents,1] = parents
	A[1:nparents,2:4] = 0
	A = as.data.frame(A)
	##########geting F1
	count_x = 1
	count_cross = 1
	B = data.frame()
	for(p1 in 1:length(P1))
	{
	    for (p2 in 1:length(P2))
	    {
	        for (i in 1:N)
	        {
	            B[count_x,1] = paste0(id,"_",count_x,"_",count_cross,"_","NA","_",i)
	            B[count_x,2] = 1
                B[count_x,3] = P1[p1]
                B[count_x,4] = P2[p2]
                B[count_x,5] = count_cross
                B[count_x,6] = NA
                B[count_x,7] = i
                count_x=count_x+1
	        }
	        count_cross = count_cross+1
	    }
	}
    C = rbind(A,as.data.frame(as.matrix(B)))
	colnames(C) = c('id','gen','p1','p2','C','S','N')
	return(C)
}

