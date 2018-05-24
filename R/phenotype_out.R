phenotype_out <- function(pheno,file,im_type)
{
    if(missing(pheno))
    {
        stop("ERROR, a phenotype object should be provided")
    }
    if(missing(file))
    {
        stop("ERROR, an output file should be provided")
    }
    if(missing(im_type))
    {im_type = "svg"}else
    {im_type = im_type}
    write.table(pheno$phenotypes,paste0(file,'_phenotypic_values.txt'),col.names=FALSE,row.names=FALSE,sep ='\t',quote = FALSE) 
    com1 = paste0(im_type,"(",file,"_phenotypic_values.",im_type,")")
    com2 = gsub("[)]","\")",gsub("[(]","(\"",com1))
    eval(parse(text = com2))
    phenotypic_values = as.numeric(as.character(pheno$phenotypes$value))
    hist(phenotypic_values,col="grey")
    dev.off()
    if("effects" %in% names(pheno))
    {
        write.table(pheno$effects,paste0(file,'_effect.txt'),col.names=FALSE,row.names=FALSE,sep ='\t',quote = FALSE)
    }else if ("training_sets" %in% names(pheno))
    {
        write.table(pheno$training_sets$train_geno,paste0(file,'_trainGeno.txt'),col.names=TRUE,row.names=TRUE,sep ='\t',quote = FALSE)
        write.table(pheno$training_sets$train_pheno,paste0(file,'_trainPheno.txt'),col.names=TRUE,row.names=TRUE,sep ='\t',quote = FALSE)
    }
}

