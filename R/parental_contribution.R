parental_contribution <- function(pop,parental_genotypes,parental_phenotypes)
{   
    founders = colnames(parental_genotypes[[1]])
    founders = data.frame(names=founders, id=seq(1,2*(length(founders)),by=2))
	rownames(founders) = founders[,1]
	col1 = colorRampPalette(c("green", "red"))(length(founders[,1])*2)
	col_founders=cbind(seq(1,length(col1)),col1)
	count = 1
	A3 = data.frame()
	for (i in 1:length(founders[,1]))
	{
		A3[count,1] = paste0(founders$name[i],'_X1')
		count = count+1
		A3[count,1] = paste0(founders$name[i],'_X2')
		count = count+1
	}
	col_founders = cbind(col_founders,A3)
	colnames(col_founders) = c('breaks','col','id')
	founders = colnames(parental_genotypes[[1]])
	founders = data.frame(names=founders, id=seq(1,2*(length(founders)),by=2))
	rownames(founders) = founders[,1]
	col1 = colorRampPalette(c("green", "red"))(length(founders[,1]))
	Col1 = cbind(id=rownames(founders),col=col1)
	haplotypes=pop
	hap_ids = names(haplotypes$haplotype)
	chr_ids = names(haplotypes$haplotype[[1]])
	scheme = pop[[1]]
	g1 = as.numeric(levels(factor(scheme$gen)))
	i = max(g1)
    A = scheme[scheme$gen==i,]
    c1 = as.numeric(levels(factor(scheme$C)))
    out = data.frame(matrix(nrow=nrow(A),ncol=length(rownames(founders))))
    colnames(out) = rownames(founders)
    rownames(out) = A$id
    for (i2 in c1)
    {
        A2 = A[A$C==i2,]
        for(x in A2$id)
        {
            B=data.frame()
            count = 1
            for (y in 1:length(chr_ids))
            {
                a = get_haplotype_matrix(haplotypes$haplotype[[x]][[y]])	
                for (z in 1:dim(a)[1])
                {
                    B[count,1] = chr_ids[y]
                    if(z == 1){B[count,2] = a[z,1]}else{B[count,2] = (a[z,1]) - (a[(z-1),1])}
                    B[count,3] = col_founders[col_founders$breaks==(a[z,2]),3]
                    B[count,4] = col_founders[col_founders$breaks==(a[z,3]),3]
                    count=count+1
                }
            }
            colnames(B) = c('chromosome_id','genetic_loci','X1','X2')
            X1 = B[,1:3]
            X2 = B[,c(1,2,4)]
            X1[,1] = paste0(X1[,1],'_X1')
            X2[,1] = paste0(X2[,1],'_X2')
            X1 = cbind(X1,group=rep('X1',nrow(X1)))
            X2 = cbind(X2,group=rep('X2',,nrow(X2)))
            colnames(X1)[3] = 'X'
            colnames(X2)[3] = 'X'
            AA = rbind(X1,X2)
            a1 = AA
            a2 = tapply(a1[,2],gsub('_.*','',a1[,3]),sum)
            out[x,which(colnames(out) %in% names(a2))] = a2
        }
    }
    out[is.na(out)] = 0
    out2 = out/apply(out,1,sum)
    values = as.matrix(t(parental_phenotypes$phenotypes$value)) 
    out3 = out2 
    for(i3 in colnames(out2))
    {
        out3[,i3] = out3[,i3] * as.numeric(as.character(values[,i3]))
    }
    out4 = apply(out3,1,sum)
    out5 = data.frame(id = names(out4) , value = out4)
    return(out5)
}

