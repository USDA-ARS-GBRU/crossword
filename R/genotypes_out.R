genotypes_out <- function(genotypes,output,level,pop,parental_genotypes,A_parent,B_parent)
{
	if(missing(level)||is.na(level)){level = NA}else
	{level = level}
	if(missing(pop) || is.na(pop)){pop = NA}else
	{pop = pop}
	do_map_format = FALSE
	if(missing(parental_genotypes) || missing(A_parent) || missing(A_parent) || A_parent =="" || B_parent == "" || is.na(A_parent) || is.na(B_parent)){do_map_format=FALSE}else
	{do_map_format = TRUE}
	if (missing(genotypes))
	{
		stop("ERROR: a genotypes object should be provided")
	}
	if(missing(output))
	{
		stop("ERROR: an output should be assinged")
	}
	geno = genotypes[[1]]
	geno2 = apply(as.data.frame(geno),1,paste0, collapse="")
	geno2b = geno2
	for (i in 1:length(geno2))
	{
		geno2b[i] = paste(sort(unlist(strsplit(geno2[i], ""))), collapse = "")
	}
	geno3 = gsub('([[:alpha:]])\\1+', '\\1', geno2b)
	for (y in 1:length(geno3))
	{
	    if(nchar(geno3[y]) ==1)
	    {
	        geno3[y] = paste0(geno3[y],geno3[y])
	    }
	}
	geno4 = sub('^(.{1})','\\1/',geno3)
	header = c(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"),colnames(geno))
	y = length(header)
	x = nrow(geno)
	geno5 = data.frame(matrix(, nrow=x, ncol=y))
	colnames(geno5) = header
	geno_c = gsub('_.*$','',rownames(geno))
	geno_p =  gsub('^.*_','',rownames(geno))
	geno5[,1] = rownames(geno)
	geno5[,2] = geno4
	geno5[,3] = geno_c
	geno5[,4] = geno_p
	geno5[,5]= '+'
	for (i in 1:ncol(geno))
	{
		geno5[,i+11] = geno[,i]
	}
	write.table(geno5,paste0(output,".hapmap"),row.names=FALSE,quote=FALSE)
	#############implementing levels
	if(!is.na(level) && !is.na(pop))
	{
        syn_geno = data.frame()
        l1 = get_level(pop[[1]],level)
        count_col = 1
        for (i in sort(unique(l1)))
        {
            s1 = geno[,which(l1 == i)] 
            for (x in 1:nrow(s1))
            {
                s2 = s1[x,]
                s3 = paste(as.matrix(s2),collapse="")
                s3_unique = rawToChar(unique(charToRaw(s3)))
                s3_unique = paste(sort(unlist(strsplit(s3_unique, ""))), collapse = "")
                if(nchar(s3_unique) != 1)
                {
                    loc1 = sub("(.).","\\1",s3_unique)
                    loc2 = sub(".(.)","\\1",s3_unique)
                    n1 = nchar(gsub(loc2,"",s3))
                    n2 = nchar(gsub(loc1,"",s3))
                    if((n1/(n1+n2)) >= 0.9)
                    {
                        s4 = paste0(loc1,loc1)
                    }else if((n2/(n1+n2)) >= 0.9)
                    {
                        s4 = paste0(loc2,loc2)
                    }else
                    {
                        s4 = paste0(loc1,loc2)
                    }
                    
                }else if(nchar(s3_unique) == 1)
                {
                    s4 = paste0(s3_unique,s3_unique)
                }
                syn_geno[x,count_col] = s4

            }
            count_col = count_col + 1
        }
        rownames(syn_geno) = rownames(geno)
        colnames(syn_geno) = sort(unique(l1))
        geno5b = cbind(geno5[,1:11],syn_geno)
        geno = syn_geno
        write.table(geno5b,paste0(output,"_synthetic.hapmap"),row.names=FALSE,quote=FALSE)
	}
	############creating biparental output for R/qtl
	if(do_map_format == TRUE)
	{
	    AA = parental_genotypes[[1]][rownames(geno),A_parent]
	    BB = parental_genotypes[[1]][rownames(geno),B_parent]
	    CC = matrix(nrow=nrow(geno),ncol=ncol(geno))
	    for(i in 1:length(AA))
	    {
	        aa = AA[i]
	        bb = BB[i]
	        uniq1 = rawToChar(unique(charToRaw(paste0(aa,collapse=""))))
            uniq2 = rawToChar(unique(charToRaw(paste0(bb,collapse=""))))
            if(nchar(uniq1) > 1 || nchar(uniq2) >1 || uniq1 == uniq2){CC[i,] = "-"}else
            {
                for(x in 1:ncol(geno))
                {
                    
                    if(geno[i,x] == aa){CC[i,x] = "A"
                    }else if(geno[i,x] == bb){CC[i,x] = "B"
                    }else if(geno[i,x] == paste0(uniq1,uniq2) || geno[i,x] == paste0(uniq2,uniq1)){CC[i,x] = "H"
                    }else{CC[i,x] = "-"}
                }
            }   
	    }
	    CC = as.data.frame(CC)
	    colnames(CC) = colnames(geno)
	    rownames(CC) = rownames(geno)
	    geno5c = cbind(geno,CC)
	    id = rownames(CC)
	    ch = gsub("_.*","",rownames(CC))
	    pos = parental_genotypes[[2]][rownames(CC),]
	    OUT = cbind(id,ch,pos,CC)
        write.table(OUT,paste0(output,"_map_format.csv"),row.names=FALSE,quote=FALSE,sep=",")
	}
	###########exporting loci information
	if(length(names(genotypes)) > 1)
	{
	    loci = genotypes[[2]]
	    loci5 = data.frame(matrix(, nrow=nrow(loci), ncol=3))
	    colnames(loci5) = c("#Chr","Pos","GenLoc")
	    loci5[,1] = geno_c
	    loci5[,2] = geno_p
	    loci5[,3] = loci[,1]
	    write.table(loci5,paste0(output,".loc"),row.names=FALSE,quote=FALSE)
	}
}

