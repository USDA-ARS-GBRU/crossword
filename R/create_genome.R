create_genome <- function(fa,genotype,output,haploid_only,max_total_size_in_gb)
{
    require("Biostrings")
    require("Rcpp")
    if(missing(haploid_only))
    {
        ho = FALSE
    }else
    {
        ho = haploid_only
    }
    if(missing(max_total_size_in_gb))
    {
        mts = 1024*1024*1024*10
    }else
    {
        mts = max_total_size_in_gb*1024*1024*1024
    }
    s = readDNAStringSet(fa)
    s2 = s
    a1 = read.table(genotype, comment.char = "", header=TRUE)
    ids = a1[,1]
    markers = a1[,12:ncol(a1)]
    row.names(markers) = ids
    expect_size = (length(colnames(markers)) * file.size(fa))
    if(mts > expect_size)
    {
        stop(paste0("ERROR: size quota is exceeded, expected size is ",round(expect_size/(mts)),"gb"))
    }
    ###cpp function to update the sequences with new SNPs
    cppFunction('List fun3(List s, List pos, List markers){
        IntegerVector Pos;
        List s2;
        s2 = s;
        std::string xx;
        CharacterVector yy;
        for(int x=0; x<s.length();x++)
        {
            Pos = pos[x];
            xx = Rcpp::as<std::string>(s2[x]);
            yy = markers[x];
            for (int y = 0;y<Pos.length();y++)
            {
                xx[Pos[y]-1] = Rcpp::as<char>(yy[y]);
            }
            s2[x] = xx;    
        }
    return(s2);}') 
    ##################################################
    dir.create(output)    
    AA = list()
    BB = list()
    CC = list()
    markers1 = list()
    markers2 = list()
    for (i2 in 1:ncol(markers))
    {
        y = colnames(markers)[i2]
        for(i in 1:length(names(s)))
        {
            AA[[names(s)[i]]] = as.character(s[[names(s)[i]]])
            BB[[names(s)[i]]] = as.matrix(markers[which(grepl(names(AA)[i],rownames(markers))),])
            CC[[names(s)[i]]] = as.numeric(cbind(pos =  gsub('^.*_','',rownames(BB[[names(s)[i]]]))))
            markers1[[names(s)[i]]] = substr(as.character(BB[[i]][,i2]),1,1)
            markers2[[names(s)[i]]] = substr(as.character(BB[[i]][,i2]),2,2) 
        }
        ss = fun3(AA,CC,markers1)
        file = file(paste0(output,"/",y,"_left_side.fa"),"w")
        for (i3 in names(ss))
        {
            writeLines(paste0(">",i3),file)
            writeLines(ss[[i3]],file)
        }
        close(file)    
        print(paste0("left_side_",colnames(markers)[i2]," was proceeded")) 
        if(ho==FALSE)
        {  
            ss = fun3(AA,CC,markers2)
            file = file(paste0(output,"/",y,"_right_side.fa"),"w")
            for (i3 in names(ss))
            {
                writeLines(paste0(">",i3),file)
                writeLines(ss[[i3]],file)
            }
            close(file)     
            print(paste0("right_side_",colnames(markers)[i2]," was proceeded"))
        }
    }
}

