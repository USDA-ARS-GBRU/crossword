## a function to estimate genomic loci from the physical loci based on sliding winows and gene frequencey
physical2genomic <- function(input,chr_stat,chr_length,window_size)
{								
	b1 = read.table(input,sep = "\t",quote = "")				#reading gff file
	b2 = b1[b1$V3=='gene',]										#extracting only gene features
	b3 = b2[order(b2$V5),]										#sorting gene ends
	b4 = b3[order(b3$V1),]										#sorting chromsomes
	rownames(b4) = gsub(';.*','',gsub('ID=','',b4$V9))			#extracting genes' IDs and naming the records
	chr_ids = levels(b4$V1)										#extracting sorted chromomes IDs
	b5 = b4[,c(1,4,5)]											#extracting dataframe of genes' IDs (record names),chromosome, gene start, gene end
	colnames(b5) = c('chr','start','end')						#naming the data frame columns
	test_cs = FALSE
	if (missing(chr_stat))										#checking if the chromosome size file is exist, otherwise the chromosome size will calculated based on the end locus (bp) of the last gene
	{
		test_cs = TRUE
	}
	else
	{
		chr_sizes = read.table(chr_stat)						#reading chromosome sizes, the file should be tab delimited containing two columns, the first column contains chromosome ID and the second contains chromsome size in bp
		colnames(chr_sizes) = c('id','size')
	}
	#checking if the chromsome length in CM was given, otherwise 100 CM will be used as default for all chromsomes
	
	if (missing(chr_length))									
	{
		cl = 100
	}
	else
	{
		if(is.numeric(chr_length))
		{
			cl = chr_length
		}
	}
	#checking if the window size was given, otherwise 100,000 bp will be used as default
	if (missing(window_size))
	{
		wz=100000
	}
	else
	{
		wz=window_size
	}
	B9 = list()													# a vector is set to hold output data structure
	if(!missing(chr_length) && is.numeric(chr_length) == FALSE && !missing(chr_stat))
	{
	    CSs = chr_sizes
	    CLs = read.table(chr_length)
	    chr_ids = intersect(as.character(CSs[,1]),as.character(CLs[,1]))
	}
	####looping through the chromosomes
	for (i in chr_ids)
	{
		b6 = b5[b5$chr==i,]										# a sub-data set containing a single chromsome data
		gene_num = dim(b6)[1]
		#reading chrosomes' sizes from the input 'chr_stat' file or estimating them from the gff file
		if(test_cs==TRUE)
		{
			chr_size=b6[dim(b6)[1],3]
		}
		else
		{
			if(i %in% chr_sizes$id==FALSE)
			{
				stop(paste0('IDs in chromosme size data does not match thoses in gff file, \ncheck chromosome \"',i),'\"')
			}
			chr_size=chr_sizes[chr_sizes$id==i,2]
		}
		
		#reading chrosomes' lengths from the input 'chr_length' file if it is not numeric.
		if(!missing(chr_length) && is.numeric(chr_length) == FALSE)
		{
			cl2 = read.table(chr_length)
			colnames(cl2) = c('id','length')
			if(i %in% cl2$id==FALSE)
			{
				stop(paste0('IDs in chromosme length data does not match thoses in gff file, \ncheck chromosome \"',i),'\"')
			}
			cl=cl2[cl2$id==i,2]
		}
		
		b8 = data.frame()										# a data frame vector to hold the tabular data
		windows = ceiling(chr_size/wz)
		x1 = 1
		y1 = 0
		y2 = 0
		#iterating through the windows
		for (i2 in 1:windows)
		{
			x2 = i2*wz
			
			b7 = subset(b6, end > x1 & end <= x2)
			b8[i2,1] = i
			b8[i2,2] = x1
			b8[i2,3] = x2
			b8[i2,4] = dim(b7)[1]/gene_num
			
			y2 = (b8[i2,4] * cl) + y1
			y1 = y2
			b8[i2,5] = y2
			
			x1 = x2
		}
		colnames(b8) = c('id','window_start','window_end','freq','length_cm')
		#data for smooth function, average of every windows, and first and last loci
		x = c(0,(b8$window_start+b8$window_end)/2,chr_size)
		y = c(0,b8$length_cm,cl)
		smo = loess(y~x)
		b9 = list(data=b8,smooth=smo,size=chr_size,length=cl)
		B9[[i]] = b9
	}
	return(B9)
}

