#the function reformates the output of cross function of simcross package
get_haplotype_matrix <- function(haplotype)
{
	a = data.frame(haplotype$mat)
	b = data.frame(haplotype$pat)	
	if ((dim(a)[1]) == 1 && (dim(b)[1]) == 1)
	{
		e=data.frame()
		e[1,1] = a[1,2]
		e[1,2] = a[1,1]
		e[1,3] = b[1,1]
	}else
	{
		a[,2] = round(a[,2],8)
		b[,2] = round(b[,2],8)
		c = as.numeric(levels(factor(c(a[,2],b[,2]))))
		c = round(c,8)
		d = data.frame()
		for (x in 1:length(c))
		{
				d[x,1] = c[x]
				
				if(c[x] %in% a[,2])
				{
					d[x,2] = a[which(c[x] == a[,2]),1]
				}
				else
				{
					d[x,2] = NA
				}
				
				if(c[x] %in% b[,2])
				{
					d[x,3] = b[which(c[x] == b[,2]),1]
				}
				else
				{
					d[x,3] = NA
				}		
		}
		e = d
		for (x in (dim(d)[1]-1):1)
		{
			if(is.na(d[x,2]))
			{
				e[x,2] = e[x+1,2]

			}
			if(is.na(d[x,3]))
			{
				e[x,3] = e[x+1,3]

			}
		}
	}
	colnames(e) = c('phy_loc','N1','N2')
	return(e)
}

