filter_parents <- function(input,output,selected_parents_file)
{
    if(missing(input))
    {
        stop("ERROR: input should be provided")
    }
    if(missing(output))
    {
        stop("ERROR: output should be provided")
    }
    if(missing(selected_parents_file))
    {
        stop("ERROR: selected_parents_file should be provided")
    }
    selected_parents = as.character(read.table(selected_parents_file)[,1])
	if (gsub('#.*','',readLines(input,n=1))=='rs')
	{
	    line = readLines(input, n = 1)
	    h2 =strsplit(line[1], ' ')[[1]]    
	    h3 = rep("NULL",length(h2))
	    h3[which(h2 %in% selected_parents)] = "character"
	    h3[1:11] = "character"
        a2 = read.table(input,comment.char="",header=TRUE,,colClasses=h3)
	    colnames(a2) = h2[h3!="NULL"]
	    write.table(a2,output,quote=FALSE,col.names=TRUE,row.names=FALSE)
	}else
	{
	    file = file(input, "r")
		while (TRUE) 
		{
			line = readLines(file, n = 1)
			h1 = strsplit(line[1], '\\t')
			if (h1[[1]][1] == "#CHROM")
			{
			  break
			}
		}
	    h2 = h1[[1]]
	    h3 = rep("NULL",length(h2))
	    h3[which(h2 %in% selected_parents)] = "character"
	    h3[1:9] = "character"
	    a2 = read.table(file,colClasses=h3)
	    colnames(a2) = h2[h3!="NULL"]
	    write.table(a2,output,quote=FALSE,col.names=TRUE,row.names=FALSE, sep = "\t")
	    close(file)    
	}
}

