add_header_button <- function(e)
{
    require(gWidgets)
    require(gWidgetstcltk)
    options(guiToolkit="tcltk")
    win1 = gwindow("crossword")
    win1 = ggroup(container = win1,use.scrollwindow=FALSE,horizontal=FALSE,width=720)
    #################input_output#####################
    ss = svalue(e)
    ss = gsub(' ','',ss)
    ss = gsub('\t*','',ss)
    ss = gsub('#.*/n','',ss)
    asd = glabel(text = "input/output",container=win1)
    font(asd) = c(size="large",weight="bold")
    g1 = ggroup(container=win1,horizontal=TRUE)
    g1a = ggroup(container=g1,horizontal=FALSE)
    g1b = ggroup(container=g1,horizontal=FALSE)    

    l1 = gedit(width = 25,text = "output_folder",container=g1a)
    t1 = gedit(width = 50,text = sub(".*output_folder=(.*)\niterations=.*","\\1",ss),container=g1b)    
    blank_t = glabel(container=g1a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b1 = gbutton(text='Browse',container=g1b,handler=function(...){
    filename <- tclvalue(tkchooseDirectory())
    if(filename!=""){svalue(t1) = paste0("\"",filename,"\"")}   
    })
    
    l20 = gedit(width = 25,text = "iterations",container=g1a)
    t20 = gedit(width = 50,text = "2",container=g1b)
       
    l2 = gedit(width = 25,text = "input_folder",container=g1a,spacing =50)
    t2 = gedit(width = 50,text = sub(".*input_folder=(.*)\nim_type=.*","\\1",ss),container=g1b,spacing=50)   
    blank_t = glabel(container=g1a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b2 = gbutton(text='Browse',container=g1b,handler=function(...){
    filename <- tclvalue(tkchooseDirectory())
    if(filename!=""){svalue(t2) = paste0("\"",filename,"\"")}   
    })
    
    l3 = gedit(width = 25,text = "im_type",container=g1a)
    rb3 = gcombobox(c("\"svg\"","\"png\"","\"pdf\"","\"jpeg\"","\"tiff\"","\"bmp\""), container=g1b,horizontal =TRUE)
    svalue(rb3) = sub(".*im_type=(.*)\nby_chromosomes=.*","\\1",ss)
    l10 = gedit(width = 25,text = "by_chromosomes",container=g1a)
    rb10 = gcombobox(c("TRUE","FALSE"),container=g1b,horizontal=TRUE)
    svalue(rb10) = sub(".*by_chromosomes=(.*)\nheterozygous=.*","\\1",ss)
    l9 = gedit(width = 25,text = "heterozygous",container=g1a)    
    rb9 = gcombobox(c("FALSE","TRUE"),container=g1b,horizontal=TRUE)
    svalue(rb9) = sub(".*heterozygous=(.*)\ngff=.*","\\1",ss)
    ###############chromosome information###########
    asd = glabel(text = "chromosome information",container=win1) 
    font(asd) = c(size="large",weight="bold")
    g2 = ggroup(container=win1,horizontal=TRUE)
    g2a = ggroup(container=g2,horizontal=FALSE)
    g2b = ggroup(container=g2,horizontal=FALSE)  
    g2c = ggroup(container=g2,horizontal=FALSE)
    l4 = gedit(width = 25,text = "gff",container=g2a)
    t4 = gedit(width = 50,text = sub(".*gff=(.*)\nchr_stat=.*","\\1",ss),container=g2b)
    blank_t = glabel(container=g2a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b3 = gbutton(text='Browse',container=g2b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t4) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l5 = gedit(width = 25,text = "chr_stat",container=g2a)
    t5 = gedit(width = 50,text = sub(".*chr_stat=(.*)\nchr_length=.*","\\1",ss),container=g2b)
    blank_t = glabel(container=g2a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b4 = gbutton(text='Browse',container=g2b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t5) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l6 = gedit(width = 25,text = "chr_length",container=g2a)
    t6 = gedit(width = 50,text = sub(".*chr_length=(.*)\nwindow_size=.*","\\1",ss),container=g2b)
    blank_t = glabel(container=g2a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b5 = gbutton(text='Browse',container=g2b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t6) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l7 = gedit(width = 25,text = "window_size",container=g2a)
    t7 = gedit(width = 50,text = sub(".*window_size=(.*)\ninput=.*","\\1",ss),container=g2b) 
    #############founder genotypes##################
    asd = glabel(text = "founder genotypes",container=win1)   
    font(asd) = c(size="large",weight="bold")
    g3 = ggroup(container=win1,horizontal=TRUE)
    g3a = ggroup(container=g3,horizontal=FALSE)
    g3b = ggroup(container=g3,horizontal=FALSE)  
    l11 = gedit(width = 25,text = "input",container=g3a)
    t11 = gedit(width = 50,text = sub(".*input=(.*)\noutcross=.*","\\1",ss),container=g3b)
    blank_t = glabel(container=g3a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b6 = gbutton(text='Browse',container=g3b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t11) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l8 = gedit(width = 25,text = "outcross",container=g3a)
    t8 = gedit(width = 50,text = sub(".*outcross=(.*)\ninput_loci=.*","\\1",ss),container=g3b)
    l11b = gedit(width = 25,text = "input_loci",container=g3a)
    t11b = gedit(width = 50,text = sub(".*input_loci=(.*)\nhomo=.*","\\1",ss),container=g3b)
    blank_t = glabel(container=g3a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b7 = gbutton(text='Browse',container=g3b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t11b) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l12 = gedit(width = 25,text = "homo",container=g3a)
    rb12 = gcombobox(c("TRUE","FALSE"),container=g3b,horizontal=TRUE)
    svalue(rb12) = sub(".*homo=(.*)\nphenotyping_method=.*","\\1",ss)
    ###############phenotyping######################
    asd = glabel(text = "phenotyping",container=win1)
    font(asd) = c(size="large",weight="bold")
    g4 = ggroup(container=win1,horizontal=TRUE)
    g4a = ggroup(container=g4,horizontal=FALSE)
    g4b = ggroup(container=g4,horizontal=FALSE)  
    l13 = gedit(width = 25,text = "phenotyping_method",container=g4a)
    rb13 <- gcombobox(c("\"QTN_random\"","\"high_low_parents\"","\"QTN_supplied\""), container=g4b,horizontal =TRUE) 
    svalue(rb13) = sub(".*phenotyping_method=(.*)\nqtn=.*","\\1",ss)
    l15 = gedit(width = 25,text = "qtn",container=g4a)
    t15 = gedit(width = 50,text = sub(".*qtn=(.*)\nh2=.*","\\1",ss),container=g4b)
    l16 = gedit(width = 25,text = "h2",container=g4a)
    t16 = gedit(width = 50,text = sub(".*h2=(.*)\nheritability_mode=.*","\\1",ss),container=g4b)
    l14 = gedit(width = 25,text = "heritability_mode",container=g4a)
    rb14 <- gcombobox(c("\"absolute\"","\"average\""), container=g4b,horizontal =TRUE)
    svalue(rb14) = sub(".*heritability_mode=(.*)\neffect_distribution=.*","\\1",ss)
    l22 = gedit(width = 25,text = "effect_distribution",container=g4a)
    rb22 = gcombobox(c("\"equal\"","\"gamma\""), container=g4b,horizontal =TRUE)
    svalue(rb22) = sub(".*effect_distribution=(.*)\nmin_qtn_freq=.*","\\1",ss)
    l23 = gedit(width = 25,text = "min_qtn_freq",container=g4a)
    t23 = gedit(width = 50,text = sub(".*min_qtn_freq=(.*)\nhighest_P=.*","\\1",ss),container=g4b)
    l14c = gedit(width = 25,text = "highest_P",container=g4a)
    t14c = gedit(width = 50,text = sub(".*highest_P=(.*)\nlowest_P=.*","\\1",ss),container=g4b)
    l14d = gedit(width = 25,text = "lowest_P",container=g4a)
    t14d = gedit(width = 50,text = sub(".*lowest_P=(.*)\nhigh_to_low_percentage=.*","\\1",ss) , container=g4b)   
    l14e = gedit(width = 25,text = "high_to_low_percentage",container=g4a)
    t14e = gedit(width = 50,text = sub(".*high_to_low_percentage=(.*)\ninput_effects=.*","\\1",ss),container=g4b)
    l19 = gedit(width = 25,text = "input_effects",container=g4a)
    t19 = gedit(width = 50,text = sub(".*input_effects=(.*)\nbiased_selection=.*","\\1",ss),container=g4b)
    blank_t = glabel(container=g4a,enabled=FALSE)
    font(blank_t) = c(size=14)
    filename=""
    dialog_b8 = gbutton(text='Browse',container=g4b,handler=function(...){
    filename <- tclvalue(tkgetOpenFile())
    if(filename!=""){svalue(t19) = paste0("\"",gsub(".*/","",filename),"\"")}   
    })
    
    l17 = gedit(width = 25,text = "biased_selection",container=g4a)
    rb17 = gcombobox(c("TRUE","FALSE"),container=g4b,horizontal=TRUE)
    svalue(rb17) = sub(".*biased_selection=(.*)\ndominant=.*","\\1",ss)
    l18 = gedit(width = 25,text = "dominant",container=g4a)
    rb18 = gcombobox(c("FALSE","TRUE"),container=g4b,horizontal=TRUE)
    svalue(rb17) = sub(".*dominant=(.*)\nI=.*","\\1",ss)
    labels = c(l1,l20,l2,l3,l10,l9,l4,l5,l6,l7,l11,l8,l11b,l12,l13,l15,l16,l14,l22,l23,l14c,l14d,l14e,l19,l17,l18)
    texts =  c(t1,t20,t2,rb3,rb10,rb9,t4,t5,t6,t7,t11,t8,t11b,rb12,rb13,t15,t16,rb14,rb22,t23,t14c,t14d,t14e,t19,rb17,rb18)
    LABELS = matrix(nrow=length(labels),ncol=1)
    TEXTS = matrix(nrow=length(labels),ncol=1)    
    g21 = ggroup(container=win1,horizontal=TRUE) 
    b1 = gbutton("apply",container=g21,handler=function(...){
    count_x = 1
    for(i in labels)
    {
        enabled(i) = FALSE
        LABELS[count_x,1] = svalue(i)
        count_x = count_x + 1
    }
    count_x = 1
    for(i in texts)
    {
        TEXTS[count_x,1] = svalue(i)
        count_x = count_x + 1
    } 
    fill_mat = paste0(LABELS,"=",TEXTS)
    svalue(e) = ""
    for (i in 1:length(labels))
    {
        add(e,fill_mat[i])
    }
    add(e,"I=\"start\"")
    if(svalue(rb13) == "\"high_low_parents\"")
    {
        if(svalue(t14c) == "" || is.na(svalue(t14c)))
        {
            stop("ERROR: you selected high_low_parents mode so you should assign the highest parent")
        }else if(svalue(t14d) == "" || is.na(svalue(t14d)))
        {
            stop("ERROR: you selected high_low_parents mode so you should assign the lowest parent")
        }
    }
    dispose(win1)
    })
    b2 = gbutton("cancel",container=g21,handler=function(...){
    dispose(win1)
    })
    editable(t19) = FALSE
    editable(t14c) = FALSE
    editable(t14d) = FALSE
    editable(t14e) = FALSE
    addHandlerClicked(rb13, handler=function(h,...) {
    if(svalue(rb13) == "\"QTN_supplied\"")
    {
        editable(t19) = TRUE
        editable(t15) = FALSE
        rb17[] = ""
        rb22[] = ""
        editable(t14c) = FALSE
        editable(t14d) = FALSE
        editable(t14e) = FALSE
    }else  if(svalue(rb13) == "\"high_low_parents\"")
    {
        editable(t19) = FALSE
        editable(t15) = TRUE
        rb17[] = c("TRUE","FALSE")
        rb22[] = c("\"equal\"","\"gamma\"")
        editable(t14c) = TRUE
        editable(t14d) = TRUE
        editable(t14e) = TRUE
        
    }else if(svalue(rb13) == "\"random\"")
    {
        editable(t19) = FALSE
        editable(t15) = TRUE
        rb17[] = c("TRUE","FALSE")
        rb22[] = c("\"equal\"","\"gamma\"")
        editable(t14c) = FALSE
        editable(t14d) = FALSE
        editable(t14e) = FALSE
    }
    })    
    source = paste0(system.file("extdata",package="crossword"),"/")
    file_1 = readLines(paste0(source,"header_info"))
    asd = gtext(container=win1,width=700)
    for (i in file_1)
    {
        if(grepl("###",i))
        {
            add(asd,gsub("### ","",i),font.attr = c(weight="bold",size="small",color = "red"))
        }else
        {
            add(asd,i,font.attr=c(size="small",color = "red"))
        }
    }
    enabled(asd) = FALSE
}

