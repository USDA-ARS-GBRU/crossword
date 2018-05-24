add_button <- function(name_function,name_output,a1,a2,e,execute_now,com,typein)
{
	require(gWidgets)
    require(gWidgetstcltk)
    options(guiToolkit="tcltk")
	w1 = gwindow(name_function,visible=FALSE,container=e)
	g0 = ggroup(horizontal=TRUE,container = w1)	
	if(name_output !="")
	{
	    l0 = glabel("output",container=g0)
	    t0 = gedit(name_output,container=g0)
	}
	if(missing(execute_now))
	{
	    en = FALSE
	}else
	{
	    en = execute_now
	}	
	if(missing(com))
	{
	    com = ""
	}	
	g1 = ggroup(container=w1,horizontal=TRUE)
    g1a = ggroup(container=g1,horizontal=F)
    g1b = ggroup(container=g1,horizontal=F)
    OUT =list()
	for (i in 1:length(a1))
	{
		eval(parse(text = paste0("l1_",i," = glabel(\"",a1[i],"\",container=g1a)")))
    if(!missing(typein))
    {
        eval(parse(text = paste0("t1_",i," = ",typein[i],"(\"",a2[i],"\",container=g1b)")))
    }else if(missing(typein))
    {
		eval(parse(text = paste0("t1_",i," = gedit(\"",a2[i],"\",container=g1b)")))
    }
    eval(parse(text =paste0("OUT[[i]] = t1_",i)))
	}
	g2 = ggroup(container=w1,horizontal=TRUE)
	out_put <- function()
	{
		if(name_output !="")
		{
		    out = paste0(svalue(t0),"=",name_function,"(")
		}else
		{
		    out = paste0(name_function,"(")
		}
		for (i in 1:length(a1))
		{
			eval(parse(text = paste0("x1 = svalue(l1_",i,")")))
			eval(parse(text = paste0("x2 = svalue(t1_",i,")")))
			out = paste0(out,x1,"=",x2,",")
		}
		out = gsub(",$",")",out)
		out = gsub(" ","",out)
		if(en==FALSE)
		{
		    add(e,out,where="at.cursor")
        }else
        {
            out = gsub("@","",out)
            eval(parse(text =out))
        }
		dispose(w1)
		return(out)		
	}
	visible(w1) <- TRUE	
	b1 = gbutton("ok",container=g2, handler = function(...) {
	out_put()} )	
	b2 = gbutton("cancel",container=g2, handler = function(...) {
	dispose(w1)} )
	if(com[[1]] != "")
	{
	    l2 = glabel(com,container=w1,wrap = FALSE,extend=TRUE)
	}
return(OUT)
}

