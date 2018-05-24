# a function is designed to visualize the relationship of physical and genomic location
# an input of this function is the output of "physical2genomic" function
# the user should specify a folder to be created and hold all plots in JPEG images
draw_PGplots <- function(input,output_folder,im_type)
{
	if (missing(im_type))
    {
        Im_type = "svg"
    }else
    {
        Im_type = im_type
    }
	dir.create(output_folder)
	for (i in 1:length(input))
	{
		if (Im_type == "svg" || Im_type == "pdf")
		{
			File = paste0(Im_type,"(\"",output_folder,'/',names(input)[i],".",Im_type,"\",20,20)")
		}else
		{
			File = paste0(Im_type,"(\"",output_folder,'/',names(input)[i],".",Im_type,"\",20,20,units=\"in\",res=300)")
		}	
		eval(parse(text = File))
		c1 = input[[i]]$data
		size = input[[i]]$size
		length = input[[i]]$length
		x = c(0,(c1$window_start+c1$window_end)/2,size)
		y = c(0,c1$length_cm,length)
		p = ggplot(data.frame(cbind(x=x/1000000,y=y)), aes(x, y)) +  geom_point() +
		xlab('Physical location in Mbp') + ylab('Genomic location in CM') + ggtitle(paste0('Chromosome ',names(input)[i])) + 
		theme(axis.text=element_text(size=40),axis.title=element_text(size=40),plot.title=element_text(size =60,face = 'bold'))
		print(p)
		dev.off()
	}
}

