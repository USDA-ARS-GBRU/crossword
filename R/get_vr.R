get_vr <- function(qtn_effect,h2,parental_genotypes,dominant,heritability_mode)
{
    if(missing(h2))
    {
        h2 = 0.5
    }else
    {
        h2 = h2
    }
    if(missing(parental_genotypes))
    {
        stop("ERROR: parental_genotypes should be provided")
    }else
    {
        parental_genotypes =parental_genotypes
    }
    if (missing(dominant))
	{
	    dominant = FALSE
	}else
	{
	    dominant = dominant
	}
	if(missing(heritability_mode))
	{
	    stop("ERROR: heritability_mode should be provided")
	}else
	{
	    heritability_mode = heritability_mode
	}
    if(is.list(qtn_effect))
    {
        Eff = qtn_effect
    }else if(is.character(qtn_effect))
    {
        ie=read.table(input_effects,header=FALSE)
        colnames(ie) = c("QTN","effect","ref")
        Eff = ie
	    Eff$effect = Eff$effect/sum(Eff$effect)
    }
    if(heritability_mode == "average")
    {
        po2 = as.data.frame(parental_genotypes$genotypes[as.character(Eff$QTN),])
        Value = get_values(po2,Eff,dominant)
        pheno = data.frame(cbind(colnames(po2),Value))
        colnames(pheno) = c('id','value')
		mat <- combn(pheno$id, 2)
		Vrs = numeric(length = ncol(mat))
		for(i in 1:ncol(mat)) {
	        l = as.character(mat[1,i])
	        h = as.character(mat[2,i])
	        geno = po2[,which(colnames(po2) %in% c(h,l))]
	        D = NA
	        for (x in 1:nrow(geno))
	        {
	            if(length(unique(as.character(geno[x,]))) > 1)
	            {
	               D = c(D,rownames(geno)[x])
	            }
	        }
	        D <- D[!is.na(D)]
	        if(length(length(D)) == 0){stop("ERROR: number of POLYMORPHIC QTNs between highest and lowest parents are less than 1, resimulate the model or increase number of QTNs")}
			Eff2 = Eff[which(as.character(Eff$QTN) %in% D),]
	        eff = as.numeric(as.character(Eff2$effect))
			#print(paste0("QTN ",qtn,"; after QTN ",length(eff)))
	        Vr = calculate_vr(eff,h2)
			#print(paste0("subVr ",Vr))
			Vrs[[i]] <- Vr
		}
		#print(Vrs)
		Vr = mean(Vrs)
		#print(paste0("meanVr ",Vr))
    }else if(heritability_mode == "absolute")
    {
        eff = as.numeric(as.character(Eff$effect))
        Vr = calculate_vr(eff,h2)
    }
    return(Vr)
}
calculate_vr <- function(eff,h2)
{
    qtn = length(eff)
    approxSize = 10000
    effMat <- matrix(eff, nrow=qtn,ncol=1)
    pMat <- matrix(sample(c(1,-1),approxSize * qtn,replace=TRUE), nrow=approxSize, ncol=qtn)
    tbv <- as.vector(pMat %*% effMat)
    vg = var(tbv)
    vr = ((vg / h2) - vg)
    return(vr)
}

