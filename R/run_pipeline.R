run_pipeline <- function(script_file)
{
    if(missing(script_file))
	{
		stop("ERROR: input script file should be provided.")
	}
	lines = readLines(script_file)
    running_script = paste0(script_file,"_running_script.R")
    write("##running R script",file=running_script,append=FALSE)
    write("library(crossword)",file=running_script,append=TRUE)
	test_g2p = FALSE
    for (i in 1:length(lines))
	{
		line=lines[i]
		line = gsub('#.*','',line)
		line = gsub(' ','',line)
		line = gsub('\t*','',line)
		if(line != '')
		{
		    if(grepl("gff=",line))
		    {
		        line = gsub("\"","",gsub("gff=","",line))
		        line = paste0(input_folder,"/",line)
		        line = paste0("gff=\"",line,"\"")
		    }
		    if(grepl("chr_stat=",line))
		    {
		        line = gsub("\"","",gsub("chr_stat=","",line))
		        line = paste0(input_folder,"/",line)
		        line = paste0("chr_stat=\"",line,"\"")
		    }
		    if(grepl("chr_length=",line))
		    {
		       line = gsub("\"","",gsub("chr_length=","",line))
		       line = paste0(input_folder,"/",line)
		       line = paste0("chr_length=\"",line,"\"")
		    }
		    if(grepl("input=",line))
		    {
		       line = gsub("\"","",gsub("input=","",line))
		       line = paste0(input_folder,"/",line)
		       line = paste0("input=\"",line,"\"")
		    }
		    if(grepl("input_loci=",line))
		    {
		       eval(parse(text = line))
		       if(!(is.na(input_loci)||input_loci==""))
		       {
		           line = gsub("\"","",gsub("input_loci=","",line))
		           line = paste0(input_folder,"/",line)
		           line = paste0("gen2phy=\"",line,"\"")
		           test_g2p = TRUE
		       }
		    }
		    if(grepl("input_pheno=",line))
		    {
		       line = gsub("\"","",gsub("input_pheno=","",line))
		       line = paste0(input_folder,"/",line)
		       line = paste0("input_pheno=\"",line,"\"")
		    }
		    if(grepl("input_effects=",line))
		    {
		       line = gsub("\"","",gsub("input_effects=","",line))
		       line = paste0(input_folder,"/",line)
		       line = paste0("input_effects=\"",line,"\"")
		    }
		    write(line,file=running_script,append=TRUE)		    
		    eval(parse(text = line))
		    if (strsplit(line,"")[[1]][1] == "I")
		    {
		        write(line,file=running_script,append=TRUE)	
		        eval(parse(text = line))
		        print("header was read")
		        write("dir.create(output_folder)",file=running_script,append=TRUE)	
		        dir.create(output_folder)
                if(test_g2p==FALSE)
                {
                    write("gen2phy = physical2genomic(gff,chr_stat,chr_length,window_size)",file=running_script,append=TRUE)
                    write("draw_PGplots(gen2phy,paste0(output_folder,\"/gen2phy_graphs\"),im_type)",file=running_script,append=TRUE)
                    write("parental_genotypes = get_parental_genotypes(input,gen2phy,homo) ",file=running_script,append=TRUE)
                    write("genotypes_out(parental_genotypes,paste0(output_folder,\"/parental_genotypes\"))",file=running_script,append=TRUE)                    
                    gen2phy = physical2genomic(gff,chr_stat,chr_length,window_size)
                    draw_PGplots(gen2phy,paste0(output_folder,"/gen2phy_graphs"),im_type)
                    parental_genotypes = get_parental_genotypes(input,gen2phy,homo) 
                    genotypes_out(parental_genotypes,paste0(output_folder,"/parental_genotypes"))
                    suppressWarnings(dput(gen2phy,file=paste0(output_folder,"/parental_genotypes_gen2phy_info"),control = "all"))
                    print("parental_genotypes were obtained")

                }else
                {
                    write("parental_genotypes = get_parental_genotypes(input,gen2phy,homo)",file=running_script,append=TRUE)                 
                    write("gen2phy = physical2genomic(gff,chr_stat,,window_size)",file=running_script,append=TRUE)
                    parental_genotypes = get_parental_genotypes(input,gen2phy,homo)
                    gen2phy = physical2genomic(gff,chr_stat,chr_length,window_size)
                }
		        break
		    }
        }
    }
    for (l in 1:iterations)
    {
        
        write("###############################################",file=running_script,append=TRUE)
        write(paste0("###iteration# ",l,"################################"),file=running_script,append=TRUE)
        write(paste0("l=",l),file=running_script,append=TRUE)
        write("dir.create(paste0(output_folder,\"/iteration_\",l))",file=running_script,append=TRUE)
        dir.create(paste0(output_folder,"/iteration_",l))
        tgv_only = FALSE
        PHENO = list()
        PHENO[["test"]] = NA
        ########getting phenotyping method
        if(phenotyping_method == "QTN_random")
        {
            write("qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution)",file=running_script,append=TRUE)
            qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution)            
        }else if(phenotyping_method == "high_low_parents")
        {
            write("qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution,highest_P=highest_P,lowest_P=lowest_P,high_to_low_percentage=high_to_low_percentage)",file=running_script,append=TRUE)
            qtn_effect = random_qtn_assign(qtn = qtn,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=dominant,effect_distribution=effect_distribution,highest_P=highest_P,lowest_P=lowest_P,high_to_low_percentage=high_to_low_percentage)   
        }else if(phenotyping_method =="QTN_supplied")
        {
            write("qtn_effect = read.table(input_effects)",file=running_script,append=TRUE)
            write("colnames(qtn_effect) = c(\"QTN\",\"ref\",\"effect\")",file=running_script,append=TRUE)
            write("rownames(qtn_effect) = qtn_effect$QTN",file=running_script,append=TRUE)            
            qtn_effect = read.table(input_effects)
            colnames(qtn_effect) = c("QTN","ref","effect")
            rownames(qtn_effect) = qtn_effect$QTN
        }
        ######getting heritability mode
        if(heritability_mode == "absolute")
        {  
            write("vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode=\"absolute\")",file=running_script,append=TRUE)            
            vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode="absolute")

        }else if(heritability_mode == "average")
        {
            write("vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode=\"average\")",file=running_script,append=TRUE)            
            vr = get_vr(qtn_effect = qtn_effect,h2 = h2 ,parental_genotypes=parental_genotypes,dominant=dominant,heritability_mode="average")
        }
        write("loci = parental_genotypes[[2]][as.character(qtn_effect$QTN),]",file=running_script,append=TRUE)        
        loci = parental_genotypes[[2]][as.character(qtn_effect$QTN),]
        write("write.table(cbind(qtn_effect,loci),paste0(output_folder,\"/iteration_\",l,\"/selected_QTNs\"),col.names=TRUE,row.names=FALSE,quote =FALSE)",file=running_script,append=TRUE)        
        write.table(cbind(qtn_effect,loci),paste0(output_folder,"/iteration_",l,"/selected_QTNs"),col.names=TRUE,row.names=FALSE,quote =FALSE)
        for (i in 1:length(lines))
	    {
	        line=lines[i]
		    line = gsub('#.*','',line)
		    line = gsub(' ','',line)
		    line = gsub('\t*','',line)
		    if(line != '')
		    {
		        if(grepl("=cross",line))
                {        
                    pop=sub("(.*)=cross.*","\\1",line)
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    Line = gsub(")$",",id=pop,chr_length=chr_length,parental_genotypes=parental_genotypes,isdh=FALSE)",line)
                    Line=gsub("=cross","=crossword::cross",Line)
                    Line=gsub("list[(]","c(",Line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: ",pop," was created"))   
                }

                if(grepl("=advance",line))
                {
                    pop=sub("(.*)=advance.*","\\1",line)
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    Line = gsub(")$",",id=pop,chr_length=chr_length)",line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: ",pop," was created"))
                }
                if(grepl("=create_families",line))
                {
                    pop=sub("(.*)=create_families.*","\\1",line)
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    Line = gsub(")$",",id=pop,chr_length=chr_length)",line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: ",pop," was created"))
                }
                if(grepl("=pick_individual",line))
                {
                    pop=sub("(.*)=pick_individual.*","\\1",line) 
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: individual ",pop," was picked"))
                }
                if(grepl("haplotypes_out",line))
                {
                    Line = gsub("output=\"",paste0("output=\"",output_folder,"/iteration_",l,"/"),line) 
                    Line = gsub(")$",paste0(",parental_genotypes=parental_genotypes)"),Line)   
                    write(Line,file=running_script,append=TRUE)                      
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: ",gsub(".*pop=(.*),output.*","\\1",Line)," was exported"))
                }
                if(grepl("=haplotypes_in",line))
                {
                    Line = gsub("=\"",paste0("=\"",output_folder,"/iteration_",l,"/"),line) 
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: individual ",sub("(.*)=haplotypes_in.*","\\",Line)," was imported"))
                }
                if(grepl("=combine_populations",line))
                {
                    Line = gsub("list[(]","c(",line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: pupulations were combined")) 
                }
                ###########################phenotyping:get_phenotype/select/reduce##################
                if(grepl("=select[(]",line) || grepl("=get_phenotype",line))
                {
                    if(grepl("tgv_only=TRUE",line))
                    {
                        tgv_only = TRUE
                    }else
                    {
                        tgv_only = FALSE
                    }
                    c1 = sub("(^.*)=.*[(].*$","\\1",line)
                    c2 = sub(".*pop=(.*)[,][Nt].*","\\1",line)
                    if(grepl("phenotype=NA",line) || grepl("phenotype=\"\"",line))
                    {
                        line = gsub(",phenotype=.*)",")",line)
                    }           
                    if(length(names(PHENO)[which(names(PHENO) %in% c2)]) == 1)
                    {
                        pheno = PHENO[[c2]]              
                    }else if (grepl("phenotype=",line))
                    {
                        com1 = sub(".*phenotype(=.*)[)]","\\1",line)
                        write(paste0("pheno", com1),file=running_script,append=TRUE)                        
                        eval(parse(text = paste0("pheno", com1)))
                        line = gsub(",phenotype=.*)",")",line)
                    }else
                    {
                        com1 = paste0("pheno = phenotype(pop=",c2,",qtn_effect=qtn_effect,tgv_only=",tgv_only,",vr=vr,parental_genotypes=parental_genotypes)")
                        write(com1,file=running_script,append=TRUE)                        
                        eval(parse(text = com1))
                        PHENO[[c2]] = pheno
                    }
                    line = gsub(",phenotype=.*)",")",line)         
                    if(grepl("=get_phenotype",line))
                    {
                        write(paste0(c1," = pheno"),file=running_script,append=TRUE)                        
                        eval(parse(text = paste0(c1," = pheno")))
                        print(paste0("Iteration_",l," #########processed: phenotyping was calculated"))
                    }else if(grepl("=select",line))
                    {
                        if(grepl("%",line))
                        {
                            percentage=TRUE
                            write("percentage=TRUE",file=running_script,append=TRUE)
                        }else
                        {
                            percentage=FALSE
                            write("percentage=FALSE",file=running_script,append=TRUE)
                        }
                        Line = gsub("%","",line)
                        sub0 = paste0("pop=",c2,",")
                        sub1 = sub(".*,(N.*)[)]","\\1",Line)
                        file_path = paste0(output_folder,"/iteration_",l,"/",c2)
                        sub2 = paste0(",pheno=pheno,percentage=percentage,im_type=im_type,file=\"",file_path,"\")")
                        com2 = paste0(c1,"=select(",sub0,sub1,sub2)
                        write(com2,file=running_script,append=TRUE)                              
                        eval(parse(text = com2))   
                        print(paste0("Iteration_",l," #########processed: a sub population was created based on selection"))
                    }
                    tgv_only <<- FALSE
                } 
                if(grepl("=dh",line))
                {
                        write(line,file=running_script,append=TRUE)                        
                        eval(parse(text = line))
                        print(paste0("Iteration_",l," #########processed: double haploid process was done"))
                }
                if(grepl("draw_haplotypes",line))
                {
                        Line = line
                        pop = sub(".*=(.*)[)]","\\1",Line)
                        write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                        ss = paste0(output_folder,"/iteration_",l,"/",pop)
                        write(paste0("ss=\"",ss,"\""),file=running_script,append=TRUE)
                        Line = gsub(")",",output_folder=ss,parental_genotypes=parental_genotypes,heterozygous=heterozygous,by_chromosomes=by_chromosomes,im_type=im_type)",Line)
                        write(Line,file=running_script,append=TRUE)                        
                        eval(parse(text = Line))
                        print(paste0("Iteration_",l," #########processed: ",pop," haplotypes were created"))              
                }    
                if(grepl("draw_population",line))
                {
                        Line = line
                        pop = sub(".*[(]pop=(.*)[)]","\\1",Line)
                        write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                        ss = paste0(output_folder,"/iteration_",l,"/",pop,"_plots")
                        write(paste0("ss=\"",ss,"\""),file=running_script,append=TRUE)                                  
                        Line = gsub(")",",file=ss,parental_genotypes=parental_genotypes,im_type=im_type)",Line)
                        write(Line,file=running_script,append=TRUE) 
                        eval(parse(text = Line))
                        print(paste0("Iteration_",l," #########processed: ",pop," graphs were created"))
                }
                if (grepl("get_genotypes",line))
                {
                    Line = gsub("get_genotypes[(]","get_genotypes(parental_genotypes=parental_genotypes,",line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: genotypes were rerieved"))
                }
                if (grepl("create_marker_set",line))
                {
                    Line = line
                    if(grepl("no_qtn=FALSE",Line))
                    {
                        Line = gsub(",no_qtn=FALSE","",Line)
                    }else if(grepl("no_qtn=TRUE",Line))
                    {
                        Line = gsub("no_qtn=TRUE","no_qtn=as.character(qtn_effect$QTN)",Line)
                    }
                    if(grepl("biased_selection=TRUE",Line))
                    {
                        Line = gsub("[)]$",",gen2phy=gen2phy,parental_genotypes=parental_genotypes)",Line)
                    }
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: the marker set was created"))
                }
                if (grepl("genomic_prediction",line))
                {
                    Line = line
                    write(Line,file=running_script,append=TRUE)
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: genomic prediction was applied"))
                }
                if (grepl("genotypes_out",line))
                {
                    Line = line
                    file_path = paste0(output_folder,"/iteration_",l)         
                    Line = gsub("output=\"",paste0("output=\"",file_path,"/"),Line)
                    Line = gsub("[)]$",",parental_genotypes=parental_genotypes)",Line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: genotypes were exported"))
                }
                if (grepl("phenotype_out",line))
                {
                    Line = line           
                    ss = sub(".*file=\"(.*)\",level.*","\\1",Line)
                    ss = paste0(output_folder,"/iteration_",l,"/",ss)               
                    ss2 = paste0(gsub(".*(,level.*)[)]","\\1",line),",N=100,method=\"top\",percentage=TRUE)")       
                    Line = gsub("file=.*[)]",paste0("file=\"",ss,"\")") ,Line)
                    Line=gsub("[)]",",im_type=im_type)",Line)
                    write(Line,file=running_script,append=TRUE)                    
                    eval(parse(text = Line))
                    file_path = paste0(output_folder,"/iteration_",l)
                    Line2 = gsub("file=.*,im_type",paste0("file=\"",ss,"\",im_type"),Line)
                    invisible(eval(parse(text = gsub("phenotype_out","select",gsub(")",ss2,Line2)))))
                    print(paste0("Iteration_",l," #########processed: the phenotypes were exported"))
                }
                if (grepl("=mas",line))
                {
                    Line = line
                    ss = sub(".*([[].*[]]).*","\\1",Line)
                    ss2 = gsub("[]]",")",gsub("[)],[(]","),c(",gsub("[[]","list(c",ss)))
                    Line = gsub("[[].*[]]",ss2,Line)
                    Line = gsub(")$",",parental_genotypes=parental_genotypes)",Line)
                    write(Line,file=running_script,append=TRUE)
                    eval(parse(text = Line))
                    print(paste0("Iteration_",l," #########processed: MAS was applied"))
                }
                if (grepl("create_population",line))
                {
                    pop = strsplit(line,'=')[[1]][1]
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    Line = gsub("@","",line)
                    Line = gsub("P=","",Line)
                    pop=sub("(^.*)=.*","\\1",Line)
                    write(paste0("pop=\"",pop,"\""),file=running_script,append=TRUE)
                    Line = gsub(")$",",chr_length=chr_length,parental_genotypes=parental_genotypes)",Line)
                    Line=gsub("list[(]","c(",Line)
                    write(Line,file=running_script,append=TRUE)
                    eval(parse(text = Line))                  
                    print(paste0("Iteration_",l," #########processed: ",pop," was created"))
                }
		    }
	    }
	    invisible(gc())
    }
}

