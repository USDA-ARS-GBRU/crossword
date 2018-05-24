phenotype <- function(pop,qtn_effect,tgv_only,vr,parental_genotypes)
{
    if(missing(pop))
    {
        stop("ERROR: pop object should be provided")
    }else
    {
        pop = pop
    }
    if(missing(qtn_effect))
    {
        stop("ERROR: qtn_effect object should be provided")
    }else
    {
        qtn_effect = qtn_effect
    }
     if(missing(vr))
    {
        stop("ERROR: vr object should be provided")
    }else
    {
        vr = vr
    }
    if(missing(tgv_only))
	{
	    tgv_only = FALSE
	}else
	{
	    tgv_only = tgv_only
	}
    ids = as.character(pop[[1]][pop[[1]]$gen==pop[[1]][nrow(pop[[1]]),2],1])
    pop2 = select_haplotype(pop,ids)
    haplo_genotypes = get_genotypes(parental_genotypes,pop)  
    po2 = as.data.frame(haplo_genotypes$genotypes[as.character(qtn_effect$QTN),])
    Value = get_values(po2,qtn_effect)
    cal_var = var(Value,na.rm = TRUE)
	cal_her = cal_var / (cal_var+vr)
    print(paste0("the calculated heritability is ",round(cal_her,3)))
	calculated_env = rnorm(ncol(po2),mean=0,sd=sqrt(vr))
    if(tgv_only == TRUE)
    {
        pheno = data.frame(cbind(colnames(po2),Value))
        colnames(pheno) = c('id','value')
    }else
    {
        pheno = data.frame(cbind(colnames(po2),Value + calculated_env, Value))	
        colnames(pheno) = c('id','value','tgv')
    }
    return(list(affected_genotypes = po2, phenotypes=pheno)) 
}

