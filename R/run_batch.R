run_batch <- function (script_file,argument_file,run)
{
    if(missing(run)){run=FALSE}else
    {run=run}
    sf = script_file
    sc = readLines(sf)
    lines = readLines(argument_file)
    list1 = NA
    list2 = NA
    for (x in 1:length(lines))
    {
        line=lines[x]
	    line = gsub('#.*','',line)
	    line = gsub(' ','',line)
	    line = gsub('\t*','',line)	
        a = line
        b = strsplit(a,"=")
        right = strsplit(b[[1]][[2]],",")[[1]]
        left = rep(b[[1]][[1]],length(right))
        if(x == 1)
        {
            list1 = left
            list2 = right
        }else
        {
         list1 = do.call(paste0,expand.grid(list1,",",left))
         list2 = do.call(paste0,expand.grid(list2,",",right))
        }
    }
    C3 = matrix(nrow=length(list2),ncol=length(strsplit(list1[1],",")[[1]]))
    out = sc[which(grepl("output_folder",sc))]
    out = gsub('#.*','',out)
	out = gsub(' ','',out)
	out = gsub('\t*','',out)	
    dir1 = gsub("\"","",gsub("output_folder=","",out))
    if(run == TRUE)
    {    
        dir.create(dir1)
    }
    for (y in 1:length(list2))
    {
        c1 = strsplit(list1[y],",")[[1]]
        c2 = strsplit(list2[y],",")[[1]]
        c3 = NA
        sc2 = sc
        for (i in 1:length(c1))
        {
            xx = sc[which(grepl(c1[i],sc))]
            xx = gsub(c1[i],c2[i],xx)
            sc2[which(grepl(c1[i],sc))] = xx
            c3 = c(c3,paste0(c1[i],"=",c2[i]))
        }
        xx = paste0("\"output_folder=combination_",y,"\"")
        sc2[which(grepl("output_folder",sc))] = gsub("^\"","", gsub("=",paste0("=\"",dir1,"/"),xx))
        c3 = c3[2:length(c3)]
        C3[y,]=c3
        file = paste0(sf,"_combination_",y)
        writeLines(sc2,paste0(sf,"_combination_",y))
        if(run == TRUE){run_pipeline(file)}
    }
    write.csv(C3,paste0(sf,"_combinations_summary.csv"))
}

