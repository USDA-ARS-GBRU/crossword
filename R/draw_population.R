draw_population <- function(pop,parental_genotypes,file,im_type)
{   
    require(mapplots)
    if (missing(im_type))
    {
        Im_type = "svg"
    }else
    {
        Im_type = im_type
    }
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
	g1 = max(g1)
	for (i in g1)
	{
	    A = scheme[scheme$gen==i,]    
        File = paste0(Im_type,"(\"",file,"_",i,".",Im_type,"\")")      
        Xlen = max(as.numeric(as.character(scheme$N)),na.rm = TRUE)       
        if(all(is.na(A$S)))
        {
            Ylen = Xlen
        }else
        {
            Ylen = max(as.numeric(as.character(A$S)),na.rm = TRUE)
        }
        Clen = max(as.numeric(as.character(A$C)),na.rm = TRUE)
        if (Im_type == "svg" || Im_type == "pdf")
		{
			 File=File
		}else
		{
			 File=gsub(")",paste0(",res=300)"),File)
		}
        eval(parse(text = File)) 
        c1 = as.numeric(levels(factor(scheme$C)))
        par(mfrow=c(Clen,1),mar=c(1,1,1,1))
        for (i2 in c1)
	    {
            plot(NA,NA, xlim=c(0,(Xlen+1)), ylim=c(0,(Ylen+1))*2,xlab=NA,ylab=NA, yaxt='n',xaxt='n')
            box(col="black",lwd=3)
            legend("bottom",legend=Col1[,1],fill=Col1[,2], horiz=TRUE)
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
                gen = A[A$id==x,2]
                C = as.numeric(as.character(A[A$id==x,5]))
                N = as.numeric(as.character(A[A$id==x,7]))
                S = as.numeric(as.character(A[A$id==x,6]))
                if(is.na(S)){S=Ylen/2}		
                add.pie(z=a2,x=N,y=S*2,labels=NA, radius=0.9,col=Col1[Col1[,1]%in%names(a2),2],edges=200)
            }
            if(!sum(is.na(A$S)) > 0)
            {
                seg = max(as.numeric(as.character(A2$S)))
                for (i3 in 2:seg)
                {
                    segments(0,((i3*2)-1),Xlen+1,((i3*2)-1),col="grey")
                }
            }
	    }
        dev.off()     
	}
}

