select <- function(pop,N,level,method,pheno,percentage,im_type,file)
{
	if(missing(pop))
	{
		stop('ERROR: haplotypes should be provided')
	}
	if(missing(N))
	{
		N=1
	}else
	{
	    N=N
	}
	if(missing(file))
	{file = "."}else
	{file = file}
	if(missing(im_type))
    {im_type = "svg"}else
    {im_type = im_type}
	if(missing(level) || level == "individual")
	{lev = 'N'}else if(level == "cross"){lev = 'C'}else if (level == "family"){lev = 'S'}else if (level == "population"){lev = "gen"}
	if(missing(method) || method == "random")
	{dec = NA}else if(method == "top"){dec = TRUE}else if(method == "bottom"){dec = FALSE}
	a = pop$scheme
	b = pop$haplotype
	for (cols in colnames(pheno$phenotypes)[2:length(colnames(pheno$phenotypes))])
	{
        pheno2 = as.numeric(as.character(pheno$phenotypes[,cols]))
        names(pheno2) = as.character(pheno$phenotypes$id)
	    if(cols == "value")
	    {
	        l1 = get_level(a,level,FALSE,TRUE)
	        l1 = l1[match(names(pheno2),a$id)]
	        a3 = data.frame()
	        for(i in sort(unique(l1)))
	        {
	            a2 = a[which(l1==i),]
	            if(!is.na(dec))
	            {
	                a2 = a2[match(names(sort(pheno2[as.character(a2$id)],decreasing = dec)),as.character(a2$id)),]
	            }
	            a3=rbind(a3,a2)
	        }
	        a = a3
            pheno3 = tapply(pheno2,l1,mean)         
            df = as.data.frame(cbind(a=pheno2,b=l1))
	        if (percentage==TRUE)
	        {
	            N = round((length(sort(unique(l1)))*N)/100)
	        }	
	        if (N > nrow(a))
	        {
		        stop('ERROR: you selected a random number equal or greater than individuals at this generation')
	        }
	        if(!is.na(dec))
	        {
	            pheno4 = sort(pheno3,decreasing = dec)
	        }else
	        {
	            pheno4 = pheno3
	        }
	        if(is.na(dec))
	        {
	            picked_pheno = names(pheno4)[sample(1:length(names(pheno4)),N)]  
	        }else
	        {
	            picked_pheno = names(pheno4)[1:N]
	        }
	        ind1 = NA
	        df2 = matrix(nrow=0,ncol=0)
	        for (i in 1:length(picked_pheno))
	        {
		        df2a = df[df$b == picked_pheno[i],]
	            if(is.na(dec))
	            {
	                df2b = df2a
	            }else
	            {
	                df2b = df2a[order(df2a$a,decreasing = dec),]
	            }
	            if(i == 1)
	            {
	                ind1 = which(df$b == picked_pheno[i])
                    df2 = df2b
	            }else
	            {
	                ind1 = c(ind1,ind1 = which(df$b == picked_pheno[i]))
                    df2 = rbind(df2,df2b)
	            }
	        }
            ids = rownames(df2)
	        out = select_haplotype(pop,ids)
	        df2 = cbind(df2,a[match(rownames(df2),as.character(a$id)),3])
	        colnames(df2) = c("value","level","mother_plant")
	        df3 = as.data.frame(pheno4[1:N])
	        colnames(df3) = "average"
	        mother_plant = NA
	        for (i2 in rownames(df3))
	        {
	            mp = as.character(df2[df2$level == i2,3])
	            if(length(sort(unique(mp))) == 1)
	            {
	                MP = sort(unique(mp))
	            }else
	            {
	                MP = NA
	            }
	            mother_plant= c(mother_plant , MP)
	        }
	        mother_plant = mother_plant[2:length(mother_plant)]
	        if(!sum(is.na(mother_plant)) > 0)
	        {
	            df3 = cbind(df3,mother_plant)
	            colnames(df3) = c("average","mother_plant")
	        }
	    }else if(cols == "tgv")
	    {
	        df2$value = pheno2[match(rownames(df2),names(pheno2))]
            df3$average = tapply(as.numeric(as.character(df2$value)),as.factor(as.character(df2$level)),mean)
	    }
	    if(cols == "value"){cols2 = "phenotypic"}else
	    {cols2 = cols}
	    write.table(df2,paste0(file,"_",level,"_",cols2,"_values.txt"),quote=FALSE)
	df3$average <- sort(df3$average)
        write.table(df3,paste0(file,"_",level,"_",cols2,"_averages.txt"),quote=FALSE)
        phenotypic_values = as.numeric(as.character(df2$value))
	}
	return(out)
}

