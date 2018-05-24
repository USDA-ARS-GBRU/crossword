simulate_fastq <- function(single_reads,input,genotypes,read_len,fold,art_binary_location)
{
    require("Biostrings")
    require("Rcpp")
    if(missing(single_reads))
    {
        single_reads = FALSE
    }else
    {
        single_reads = single_reads
    }
    if(missing(input))
    {
        stop("ERROR: input file or folder should be provided")
    }
    if(missing(read_len))
    {
        read_len = 250
    }
    if(missing(fold))
    {
        fold = 30
    }
    if(missing(art_binary_location) || is.na(art_binary_location) || art_binary_location == "")
    {
        abl = paste0(system.file("art",package="crossword"),"/",Sys.info()["sysname"])
    }else
    {
        abl = art_binary_location
    }
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
    if(dir.exists(input))
    {
        files = list.files(input)
        folds2 = fold/(length(files))
        for (i in files)
        {
            fa = paste0(input,"/",i)
            com1 = paste0(abl,"/art_illumina -i ",fa," -l ",read_len," -f ",folds2," -o ",paste0(input,"/"),"_",i)
            if(single_reads == FALSE)
            {
                com1 = gsub("$",paste0(" -m ",round(read_len *2 /1.5,0)," -s 10"),com1)
            }
            system(com1)      
        }
    }else if (is.character(input) && missing(genotypes))
    {
        fa = input
        com1 = paste0(abl, "/art_illumina -i ",fa," -l ",read_len," -f ",fold," -o ",input)
        if(single_reads == FALSE)
        {
            com1 = gsub("$",paste0(" -m ",round(read_len *2 /1.5,0)," -s 10"),com1)
        }
        system(com1)
    }else if (is.character(input) && !missing(genotypes))
    {   
        fa = input
        s = readDNAStringSet(fa)
        s2 = s
        a1 = read.table(genotypes, comment.char = "", header=TRUE)
        ids = a1[,1]
        markers = a1[,12:ncol(a1)]
        folds2 = fold / (ncol(markers)*2)
        Input = paste0(genotypes,"_simulated_reads")
        dir.create(Input)
        row.names(markers) = ids     
        AA = list()
        BB = list()
        CC = list()
        markers1 = list()
        markers2 = list()
        for (i2 in 1:ncol(markers))
        {
            for(i in 1:length(names(s)))
            {
                AA[[names(s)[i]]] = as.character(s[[names(s)[i]]])
                BB[[names(s)[i]]] = as.matrix(markers[which(grepl(names(AA)[i],rownames(markers))),])
                CC[[names(s)[i]]] = as.numeric(cbind(pos =  gsub('^.*_','',rownames(BB[[names(s)[i]]]))))
                markers1[[names(s)[i]]] = substr(as.character(BB[[i]][,i2]),1,1)
                markers2[[names(s)[i]]] = substr(as.character(BB[[i]][,i2]),2,2) 
            }
            ss = fun3(AA,CC,markers1)
            file = file("temp_532907", "w")
            for (i3 in names(ss))
            {
                writeLines(paste0(">",i3),file)
                writeLines(ss[[i3]],file)
            }
            close(file)          
            ss = fun3(AA,CC,markers2)
            file = file("temp_532908", "w")
            for (i3 in names(ss))
            {
                writeLines(paste0(">",i3),file)
                writeLines(ss[[i3]],file)
            }    
            close(file)
            
            invisible(file.append("temp_532907","temp_532908"))
            com1 = paste0(abl,"/art_illumina -i ","temp_532907"," -l ",read_len," -f ",folds2," -o ",paste0(Input,"/"),"_",colnames(markers)[i2],"_")
            if(single_reads == FALSE)
            {
                com1 = gsub("$",paste0(" -m ",round(read_len *2 /1.5,0)," -s 10"),com1)
            }
            system(com1)
            
            file.remove("temp_532907")
            file.remove("temp_532908")
            
            print(paste0(colnames(markers)[i2],"was proceeded"))
        }
    }
}

