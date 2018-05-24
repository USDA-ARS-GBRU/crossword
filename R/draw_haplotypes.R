draw_haplotypes <- function(haplotypes,output_folder,parental_genotypes,heterozygous,by_chromosomes,im_type)
{
	require('gridExtra')
	require('ggplot2')
	if (missing(im_type))
    {
        Im_type = "svg"
    }else
    {
        Im_type = im_type
    }
	if(missing(haplotypes))
	{
		stop('ERROR: input haplotypes should be provided')
	}
	if(missing(output_folder))
	{
		stop('ERROR: output_folder file prefex should be provided')
	}
	if(missing(output_folder))
	{
		stop('ERROR: founder list should be provided')
	}
	het = 0
	if(missing(heterozygous))
	{
		het = 0
	}else if(heterozygous == TRUE)
	{
		het = 1
	}
	sty = 0
	if(missing(by_chromosomes))
	{
		sty = 0
	}else if(by_chromosomes == TRUE)
	{
		sty = 1
	}
	dir.create(output_folder)
	founders = colnames(parental_genotypes[[1]])
	founders = data.frame(names=founders, id=seq(1,2*(length(founders)),by=2))
	rownames(founders) = founders[,1]
	col1 = colorRampPalette(c("green", "red"))((length(founders[,1]))*2)
	if(het == 0){for (i in seq(2,length(col1),by = 2)){col1[i] = col1[i-1]}}
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
	sch = haplotypes[[1]]
	hap_ids = sch[sch$gen == sch[nrow(sch),2],1]
    haplotypes = select_haplotype(haplotypes,hap_ids)	
	chr_ids = names(haplotypes$haplotype[[1]])
	if (sty == 1)
	{
		for (x in 1:length(hap_ids))
		{
			A=data.frame()
			count = 1
			for (y in 1:length(chr_ids))
			{
				a = get_haplotype_matrix(haplotypes$haplotype[[x]][[y]])	
				for (z in 1:dim(a)[1])
				{
					#hap_ids[x]
					A[count,1] = chr_ids[y]
					if(z == 1){A[count,2] = a[z,1]}else{A[count,2] = (a[z,1]) - (a[(z-1),1])}
					A[count,3] = col_founders[col_founders$breaks==(a[z,2]),3]
					A[count,4] = col_founders[col_founders$breaks==(a[z,3]),3]
					count=count+1
				}
			}
			colnames(A) = c('chromosome_id','genetic_loci','X1','X2')
			X1 = A[,1:3]
			X2 = A[,c(1,2,4)]
			X1[,1] = paste0(X1[,1],'_X1')
			X2[,1] = paste0(X2[,1],'_X2')
			X1 = cbind(X1,group=rep('X1',nrow(X1)))
			X2 = cbind(X2,group=rep('X2',,nrow(X2)))
			colnames(X1)[3] = 'X'
			colnames(X2)[3] = 'X'
			AA = rbind(X1,X2)
			PR = factor(AA$X) #parental_recombination
			col2 = t(as.matrix(col_founders[col_founders$id %in%  levels(factor(AA$X)),2:3]))
			colnames(col2) = col2[2,]
			col2 = col2[1,]
			if (Im_type == "svg" || Im_type == "pdf")
			{
				File = paste0(Im_type,"(\"",output_folder,'/',hap_ids[x],".",Im_type,"\",y*1.5,10)")
			}else
			{
				File = paste0(Im_type,"(\"",output_folder,'/',hap_ids[x],".",Im_type,"\",y*1.5,10,units=\"in\",res=300)")
			}
			eval(parse(text = File))	
			g1 = ggplot(data  = AA, aes(x = chromosome_id, y = genetic_loci,group = factor(group), fill = PR)) + geom_bar(stat = "identity") +  
				xlab('Chromsome_IDs') + ylab ('Genetic_loci (CM)') + ggtitle(hap_ids[x]) + theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5,size=30),axis.title=element_text(size=40),plot.title=element_text(size =60,face = 'bold')) +
				ggtitle(hap_ids[x]) + scale_fill_manual(values = col2)
			print(g1)
			dev.off()	
		}
	}else if (sty == 0)
	{
		for (y in 1:length(chr_ids))
		{
			A=data.frame()
			count = 1
			for (x in 1:length(hap_ids))
			{
				a = get_haplotype_matrix(haplotypes$haplotype[[x]][[y]])	
				for (z in 1:dim(a)[1])
				{
					A[count,1] = hap_ids[x]
					if(z == 1){A[count,2] = a[z,1]}else{A[count,2] = (a[z,1]) - (a[(z-1),1])}
					A[count,3] = col_founders[col_founders$breaks==(a[z,2]),3]
					A[count,4] = col_founders[col_founders$breaks==(a[z,3]),3]
					count=count+1
				}
			}
			colnames(A) = c('individual_id','genetic_loci','X1','X2')
			X1 = A[,1:3]
			X2 = A[,c(1,2,4)]			
			X1 = cbind(X1,group=rep('X1',nrow(X1)))
			X2 = cbind(X2,group=rep('X2',,nrow(X2)))
			colnames(X1)[3] = 'X'
			colnames(X2)[3] = 'X'
			AA = rbind(X1,X2)
			PR = factor(AA$X)
			col2 = t(as.matrix(col_founders[col_founders$id %in%  levels(factor(AA$X)),2:3]))
			colnames(col2) = col2[2,]
			col2 = col2[1,]
			File = paste0(Im_type,"(\"",output_folder,'/',chr_ids[y],".",Im_type,"\",x*1.5,10)")
			eval(parse(text = File))
			g1 = ggplot(data  = AA, aes(x = individual_id, y = genetic_loci,group = factor(group), fill = PR)) + geom_bar(stat = "identity") +  
			ggtitle(chr_ids[y]) + 
			xlab('Individual_IDs') + ylab ('Genetic_loci (CM)') + 
			theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5,size=30),axis.title=element_text(size=40),plot.title=element_text(size =60,face = 'bold') ) +
			scale_fill_manual(values = col2) 
			print(g1)
			dev.off()	
		}
	}
}

