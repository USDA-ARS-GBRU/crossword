crossword_gui <- function()
{
    require(gWidgets)
    require(gWidgetstcltk)
    options(guiToolkit="tcltk")
    win0 = gwindow("crossword",height=1500,width=2500)
    t0 = gtext(container=win0,wrap=TRUE,extend=TRUE)
    source = paste0(system.file("extdata",package="crossword"),"/")
    script_file =  paste0(source,"crossword_script_input_file_gui.script")
    lines = readLines(script_file)
    file_2 = readLines(paste0(source,"gui_comments"))
    svalue(t0) <- lines
    add(t0,"")
    list1 = list()
    ##header
    list1$header$edit_header$handler = function(h,...)
        {add_header_button(t0)}
    ###input_output
    list1$input_output$haplotypes_in$handler =  function(h,...)
        {add_button("haplotypes_in","pop5",c("input"),c("\\\"haplo1_out\\\""),t0,FALSE,com = comments(file_2,"haplotypes_out"))}
    list1$input_output$haplotypes_out$handler = function(h,...)
        {add_button("haplotypes_out","",c("pop","output"),c("pop2","\\\"haplo1_out\\\""),t0,FALSE,com = comments(file_2,"haplotypes_in"))} 
    list1$input_output$draw_haplotypes$handler = function(h,...)
        {add_button("draw_haplotypes","",c("haplotypes"),c("pop2"),t0,com = comments(file_2,"draw_haplotypes"))}
    list1$input_output$draw_population$handler = function(h,...)
        {add_button("draw_population","",c("pop"),c("pop3"),t0,FALSE,com = comments(file_2,"draw_population"))}
    list1$input_output$phenotype_out$handler = function(h,...)
        {po2 = add_button("phenotype_out","",c("pheno","file","level","pop"),c("pheno3", "\\\"pheno3_predicted\\\"","\\\"family\\\"","pop3"),t0,FALSE,com = comments(file_2,"phenotype_out"),c("gedit","gedit","gcombobox","gedit"))
        po3 = po2[[3]]
        po3[] = c("NA","\"individual\"","\"family\"","\"cross\"","\"population\"")
        }
    list1$input_output$genotypes_out$handler = function(h,...)
        {go2 = add_button("genotypes_out","",c("genotypes","output","level","pop","A_parent","B_parent"),c("ms3","\\\"pop3_marker_set\\\"","\\\"family\\\"","pop3","NA","NA"),t0,FALSE,com = comments(file_2,"genotype_out"),c("gedit","gedit","gcombobox","gedit","gedit","gedit"))
        go3 = go2[[3]]
        go3[] = c("NA","\"individual\"","\"family\"","\"cross\"","\"population\"")
        }    
    ###population_manipulation
    list1$population_manipulation$cross$handler = function(h,...)
        {add_button("cross","pop1",c("P1","P2","N"),c("\\\"grg\\\"","list(\\\"tr\\\",\\\"tg\\\")","10"),t0,FALSE,com = comments(file_2,"cross"))}
    list1$population_manipulation$advance$handler = function(h,...)
        {adv2 = add_button("advance","pop2",c("pop","F","level","clevel"),c("pop1","5","0","0"),t0,FALSE,com = comments(file_2,"advance"),c("gedit","gedit","gcombobox","gcombobox"))
        adv3 = adv2[[3]]
        adv3[] = c("\"individual\"","\"family\"","\"cross\"","\"population\"")
        adv4 = adv2[[4]]
        adv4[] = c("\"individual\"","\"family\"","\"cross\"","\"population\"")
        }
    list1$population_manipulation$create_families$handler = function(h,...)
        {cf2 = add_button("create_families","pop3",c("pop","S","is_clone"),c("pop2","5","FALSE"),t0,FALSE,com = comments(file_2,"create_families"),c("gedit","gedit","gcombobox"))
        cf3 = cf2[[3]]
        cf3[] = c("FALSE","TRUE")}   
    list1$population_manipulation$select$handler = function(h,...)
        {select2 = add_button("select","pop6",c("pop","N","level","method","phenotype"),c("pop3","5%","0","\\\"top\\\"","NA"),t0,FALSE,com = comments(file_2,"select"),c("gedit","gedit","gcombobox","gcombobox","gedit"))
        select3 = select2[[3]]
        select3[] = c("\"individual\"","\"family\"","\"cross\"","\"population\"")
        select4 = select2[[4]]
        select4[] = c("\"top\"","\"bottom\"","\"random\"") } 
    list1$population_manipulation$dh$handler = function(h,...)
        {add_button("dh","pop10",c("pop"),c("pop1"),t0,FALSE,com = comments(file_2,"dh"))}
    list1$population_manipulation$mas$handler = function(h,...)
        {mas2 = add_button("mas","pop3_b",c("pop","marker","level"),c("pop3","[(\\\"Aradu.A01_8315233\\\",\\\"A\\\",0.5),(\\\"Araip.B10_125281670\\\",\\\"C\\\",0.5)]","1"),t0,FALSE,com = comments(file_2,"mas"),c("gedit","gedit","gcombobox"))
        mas3 = mas2[[3]]
        mas3[] = c("\"individual\"","\"family\"","\"cross\"","\"population\"")
        }
    list1$population_manipulation$create_population$handler = function(h,...)
        {add_button("create_population","pop10",c("P"),c("list(\\\"tg\\\",\\\"tr\\\",\\\"grg\\\")"),t0,FALSE,com = comments(file_2,"create_population"))}
    list1$population_manipulation$pick_individual$handler = function(h,...)
        {add_button("pick_individual","pop11",c("pop","cross","family","individual"),c("pop10","1","2","15"),t0,FALSE,com = comments(file_2,"pick_individual"))}  
    ###helper
    list1$helper$combine_populations$handler = function(h,...)
        {add_button("combine_populations","pop9","pops","list(pop1,pop2)",t0,FALSE,com = comments(file_2,"combine_populations"))}
    list1$helper$get_phenotype$handler = function(h,...)
        {gp2 = add_button("get_phenotype","pheno2",c("pop","tgv_only"),c("pop2","FALSE"),t0,FALSE,com = comments(file_2,"get_phenotype"),c("gedit","gcombobox"))
        gp3 = gp2[[2]]
        gp3[] = c("FALSE","TRUE")
        }
    list1$helper$genomic_prediction$handler = function(h,...)
        {gs2 =add_button("genomic_prediction","pheno3",c("train_geno","train_pheno","predict_geno_pop","method","level","train_pop"),c("ms2","pheno2","ms3","\\\"rrBLUP\\\"","\\\"cross\\\"","pop2"),t0,FALSE,com = comments(file_2,"genomic_prediction"),c("gedit","gedit","gedit","gcombobox","gcombobox","gedit"))
        gs3 = gs2[[4]]
        gs3[] = c("\"rrBLUP\"","\"parental_contribution\"")
        gs4 = gs2[[5]]
        gs4[] = c("\"individual\"","\"family\"","\"cross\"","\"population\"")
        }
    list1$helper$get_genotypes$handler = function(h,...)
        {add_button("get_genotypes","geno2",c("pop","pre_selected_markers"),c("pop3","NA"),t0,FALSE,com = comments(file_2,"get_genotypes"))} 
    list1$helper$create_marker_set$handler = function(h,...)
        {cms2 = add_button("create_marker_set","ms2",c("genotypes","N","MAF","no_qtn","biased_selection"),c("geno2","200","0.13","","FALSE"),t0,FALSE,com = comments(file_2,"create_marker_set"),c("gedit","gedit","gedit","gcombobox","gcombobox"))
        cms3 = cms2[[4]]
        cms3[] = c("FALSE","TRUE")
        cms4 = cms2[[5]]
        cms4[] = c("FALSE","TRUE")
        }    
    ###running    
    list1$applying$apply$handler =  function(h,...) 
        {
            file = tclvalue(tkgetSaveFile(initialfile="./temp_script_file"))
            writeLines(svalue(t0),file)
            if(is.na(file) || file == "")
            {
                gmessage(message="ERROR: you have to assign an input file",title="ERROR",icon="error")
                stop("ERROR: you have to assign an input file")
            }else
            {
                ptm = proc.time()
                run_pipeline(file)
                st1 = round((proc.time() - ptm)[3]/60,2)
                gmessage(message=paste0("the processing time was ",st1," minutes"),title="processing time",icon="info")
            }
        }
    list1$applying$open$handler =  function(h,...) 
        {
            file = tclvalue(tkgetOpenFile(initialfile="./temp_script_file"))
            if(is.na(file) || file == "")
            {
                gmessage("message=ERROR: you have to assign an input file",title="ERROR",icon="error")
                stop("ERROR: you have to assign an input file")
            }else
            {
                svalue(t0) = readLines(file)
            }
        }
    list1$applying$cancel$handler =  function(h,...)    
        {dispose(win0)}       
    ###auxiliary_programs
    list1$auxiliary_programs$count_chr_length$handler =  function(h,...) 
        {add_button("count_chr_length","",c("fa","output"),c("\\\"input.fasta\\\"","\\\"chromosome_length\\\""),t0,TRUE,com = comments(file_2,"count_chr_legnth"))}
    list1$auxiliary_programs$vcf2hapmap$handler =  function(h,...) 
        {add_button("vcf2hapmap","",c("input","output"),c("\\\"input.vcf\\\"","\\\"output.hapmap\\\""),t0,TRUE,com = comments(file_2,"vcf2hapmap"))}
    list1$auxiliary_programs$create_genome$handler =  function(h,...) 
        {cg2 = add_button("create_genome","",c("fa","genotype","output","haploid_only","max_total_size_in_gb"),c("\\\"input.fasta\\\"","\\\"geno3.hapmap\\\"","\\\"output\\\"","FALSE","20"),t0,TRUE,com = comments(file_2,"create_genome"),c("gedit","gedit","gedit","gcombobox","gedit"))
        cg3 = cg2[[4]]
        cg3[] = c("FALSE","TRUE")
        }
    list1$auxiliary_programs$simulate_fastq$handler =  function(h,...) 
        {sim2 = add_button("simulate_fastq","",c("single_reads","input","genotypes","read_len","fold","art_binary_location"),c("FALSE","\\\"input.fasta\\\"","\\\"geno3.vcf\\\"","250","10",""),t0,TRUE,com = comments(file_2,"simulate_fastq"),c("gcombobox","gedit","gedit","gedit","gedit","gedit"))
        sim3 = sim2[[1]]
        sim3[] = c("FALSE","TRUE")
        }
    list1$auxiliary_programs$filter_parents$handler =  function(h,...) 
        {add_button("filter_parents","",c("input","output","selected_parents_file"),c("\\\"peanut.vcf\\\"","\\\"filtered_peanut.vcf\\\"","\\\"selected_parents_file\\\""),t0,TRUE,com = comments(file_2,"filter_parents"))}  
    ######################    
    mb <- gmenu(list1, container=win0)
    advance1 = paste0(system.file("icons",package="crossword"),"/advance.gif")
    combine_populations1 = paste0(system.file("icons",package="crossword"),"/combine_populations.gif")
    convert1 = paste0(system.file("icons",package="crossword"),"/convert.gif")
    count1 = paste0(system.file("icons",package="crossword"),"/count.gif")
    create_families1 = paste0(system.file("icons",package="crossword"),"/create_families.gif")
    create_genome1 = paste0(system.file("icons",package="crossword"),"/create_genome.gif")
    create_populations1 = paste0(system.file("icons",package="crossword"),"/create_populations.gif")
    cross1 = paste0(system.file("icons",package="crossword"),"/cross.gif")
    dh1 = paste0(system.file("icons",package="crossword"),"/dh.gif")
    draw1 = paste0(system.file("icons",package="crossword"),"/draw.gif")
    edit_header1 = paste0(system.file("icons",package="crossword"),"/edit_header.gif")
    filter1 = paste0(system.file("icons",package="crossword"),"/filter.gif")
    genomic_prediction1 = paste0(system.file("icons",package="crossword"),"/genomic_prediction.gif")
    get1 = paste0(system.file("icons",package="crossword"),"/get.gif")
    in1 = paste0(system.file("icons",package="crossword"),"/in.gif")
    mas1 = paste0(system.file("icons",package="crossword"),"/mas.gif")
    out1 = paste0(system.file("icons",package="crossword"),"/out.gif")
    select1 = paste0(system.file("icons",package="crossword"),"/select.gif")
    simulate1 = paste0(system.file("icons",package="crossword"),"/simulate.gif")
    ####################
    a0 = gaction(label = "external_files_phenotyping >>>>> ", handler=function(h,...) NA)
    a1 = gaction(label = names(list1[[1]])[1], icon=edit_header1,  handler=list1$header[[1]]$handler)
    tl1 = list(header=a0,edit_header=a1)
    b0 = gaction(label = "input_output >>>>> ", handler=function(h,...) NA)
    b1 = gaction(label = names(list1[[2]])[1], icon=in1,  handler=list1$input_output[[1]]$handler)
    b2 = gaction(label = names(list1[[2]])[2], icon=out1,  handler=list1$input_output[[2]]$handler)
    b3 = gaction(label = names(list1[[2]])[3], icon=draw1,  handler=list1$input_output[[3]]$handler)
    b4 = gaction(label = names(list1[[2]])[4], icon=draw1,  handler=list1$input_outpu[[4]]$handler)
    b5 = gaction(label = names(list1[[2]])[5], icon=out1,  handler=list1$input_output[[5]]$handler)
    b6 = gaction(label = names(list1[[2]])[6], icon=out1,  handler=list1$input_output[[6]]$handler)
    tl2 <- list(input_output = b0,haplotypes_in=b1,haplotypes_out=b2,draw_haplotypes=b3,draw_population=b4,phenotype_out=b5,genotypes_out=b6)
    c0 = gaction(label = "population_manipulation >>>>> ", handler=function(h,...) NA)
    c1 = gaction(label = names(list1[[3]])[1], icon=cross1,  handler=list1$population_manipulation[[1]]$handler)
    c2 = gaction(label = names(list1[[3]])[2], icon=advance1,  handler=list1$population_manipulation[[2]]$handler)
    c3 = gaction(label = names(list1[[3]])[3], icon=create_families1,  handler=list1$population_manipulation[[3]]$handler)
    c4 = gaction(label = names(list1[[3]])[4], icon=select1,  handler=list1$population_manipulation[[4]]$handler)
    c5 = gaction(label = names(list1[[3]])[5], icon=dh1,  handler=list1$population_manipulation[[5]]$handler)
    c6 = gaction(label = names(list1[[3]])[6], icon=mas1,  handler=list1$population_manipulation[[6]]$handler)
    c7 = gaction(label = names(list1[[3]])[7], icon=create_populations1,  handler=list1$population_manipulation[[7]]$handler)
    c8 = gaction(label = names(list1[[3]])[8], icon=select1,  handler=list1$population_manipulation[[8]]$handler)
    tl3 =  list(population_manipulation = c0,cross=c1,advance=c2,create_families=c3,select=c4,dh=c5, mas=c6, create_population=c7,select_individuals=c8)
    d0 = gaction(label = "helper >>>>> ", handler=function(h,...) NA)
    d1 = gaction(label = names(list1[[4]])[1], icon=combine_populations1,  handler=list1$helper[[1]]$handler)
    d2 = gaction(label = names(list1[[4]])[2], icon=get1,  handler=list1$helper[[2]]$handler)
    d3 = gaction(label = names(list1[[4]])[3], icon=genomic_prediction1,  handler=list1$helper[[3]]$handler)
    d4 = gaction(label = names(list1[[4]])[4], icon=get1,  handler=list1$helper[[4]]$handler)
    d5 = gaction(label = names(list1[[4]])[5], icon=mas1,  handler=list1$helper[[5]]$handler)
    tl4 = list(helper=d0,combine_populations=d1,c7,get_phenotype=d2,genomic_prediction = d3,get_genotypes=d4,create_marker_set=d5)
    e0 = gaction(label = "manage >>>>> ", handler=function(h,...) NA)
    e1 = gaction(label = names(list1[[5]])[1], icon="ok",  handler=list1$applying[[1]]$handler)
    e2 = gaction(label = names(list1[[5]])[2], icon="open",  handler=list1$applying[[2]]$handler)
    e3 = gaction(label = names(list1[[5]])[3], icon="cancel",  handler=list1$applying[[3]]$handler)
    tl5 = list(applying=e0,apply=e1,open=e2,cancel=e3)
    f0 = gaction(label = "auxiliary_programs >>>>> ", handler=function(h,...) NA)
    f1 = gaction(label = names(list1[[6]])[1], icon=count1,  handler=list1$auxiliary_programs[[1]]$handler)
    f2 = gaction(label = names(list1[[6]])[2], icon=convert1,  handler=list1$auxiliary_programs[[2]]$handler)
    f3 = gaction(label = names(list1[[6]])[3], icon=create_genome1,  handler=list1$auxiliary_programs[[3]]$handler)
    f4 = gaction(label = names(list1[[6]])[4], icon=simulate1,  handler=list1$auxiliary_programs[[4]]$handler)
    f5 = gaction(label = names(list1[[6]])[5], icon=filter1,  handler=list1$auxiliary_programs[[5]]$handler)
    tl6 = list(auxiliary_programs = f0,count_chr_length=f1,vcf2hapmap=f2,create_genome=f3,simulate_fastq=f4,filter_parents=f5)
    ##############
    gtoolbar(tl1, container = win0)
    gtoolbar(tl2, container = win0)
    gtoolbar(tl3, container = win0)
    gtoolbar(tl4, container = win0)
    gtoolbar(tl5, container = win0)
    gtoolbar(tl6, container = win0)
    enabled(tl1[[1]]) = FALSE
    enabled(tl2[[1]]) = FALSE
    enabled(tl3[[1]]) = FALSE
    enabled(tl4[[1]]) = FALSE
    enabled(tl5[[1]]) = FALSE
    enabled(tl6[[1]]) = FALSE
}
######function to create comments
comments <- function(file,function_name)
{
    b = file[(which(file == paste0("### ",function_name))+1):length(file)]
    c = b[1:(which(b == b[grep("###.*",b)][1])-1)]
    return(c)
}

