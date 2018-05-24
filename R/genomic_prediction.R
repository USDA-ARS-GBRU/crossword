genomic_prediction <- function(train_geno,train_pheno,predict_geno_pop,method,level,train_pop)
{
    require("rrBLUP")
    if(missing(method))
    {
        M = "rrBLUP"
    }else if (method == "rrBLUP")
    {
        M = "rrBLUP"
    }else if (method == "parental_contribution")
    {
        M = "parental_contribution"
    }else
    {
        stop("ERROR: method is not available")
    }
    if(missing(level) || level == "" || is.na(level))
    {
        level = NA
    }else if (level == "individual")
    {
        level = "individual"
    }else
    {
        level = level
    }
    if(missing(train_pop) || train_pop == "" || is.na(train_pop))
    {
        train_pop = NA
    }else
    {
        train_pop = train_pop
    }
    TrainGeno = NA
    TrainPheno = NA
    rrplub_geno <- function(genotype)
    {
        A = data.frame()
        a1 = genotype
        for (x in 1:nrow(a1))
        {
            uniq1 = rawToChar(unique(charToRaw(paste0(as.matrix(a1[x,]),collapse=""))))
            uniq1 = paste(sort(unlist(strsplit(uniq1, ""))), collapse = "")
            snp1 = strsplit(uniq1,"")[[1]][1]
            snp2 = strsplit(uniq1,"")[[1]][2]
            for (y in 1:ncol(a1))
            {
                if(nchar(uniq1) == 1)
                {
                    A[x,y] = 0
                }else
                {
                    if(a1[x,y] == paste0(snp1,snp1))
                    {
                        A[x,y] = 1
                    }else if(a1[x,y] == paste0(snp2,snp2))
                    {
                        A[x,y] = -1
                    }else if(a1[x,y] == paste0(snp1,snp2))
                    {
                        A[x,y] = 0
                    }else if(a1[x,y] == paste0(snp2,snp1))
                    {
                        A[x,y] = 0
                    }
                }
            }
        }
        colnames(A) = colnames(a1)
        rownames(A) = rownames(a1)
        A = t(A)
        return(A)
    }  
    predict_pheno = NA
    a2 = rrplub_geno(predict_geno_pop[[1]])
    if (M == "rrBLUP")
    {
        if((is.na(level) && is.na(train_pop)) || level == "individual")
        {
            a1 = rrplub_geno(train_geno[[1]])
            b1 = as.double(as.character(train_pheno$phenotypes$value))
            ans = kinship.BLUP(y=as.vector(b1),G.train=a1,G.pred=a2,K.method="GAUSS")
            TrainGeno = train_geno[[1]]
            TrainPheno = b1

        }else if(level == "cross" || level == "family" || level == "population")
        {
            l1 = get_level(train_pop$scheme,level)
            train_phen_values = as.double(as.character(train_pheno$phenotypes$value))
            pheno = tapply(train_phen_values,l1,mean)
            m2 = matrix(ncol=2,nrow = length(pheno))
            m2[,1] = as.matrix(pheno)
            syn_geno = data.frame()
            count_col = 1
            for (i in sort(unique(l1)))
            {
                s1 = as.matrix(train_geno[[1]])[,which(l1 == i)] 
                for (x in 1:nrow(s1))
                {
                    s2 = s1[x,]
                    s3 = paste(s2,collapse="")
                    s3_unique = rawToChar(unique(charToRaw(s3)))
                    if(nchar(s3_unique) != 1)
                    {
                        loc1 = sub("(.).","\\1",s3_unique)
                        loc2 = sub(".(.)","\\1",s3_unique)
                        n1 = nchar(gsub(loc2,"",s3))
                        n2 = nchar(gsub(loc1,"",s3))

                        if((n1/(n1+n2)) >= 0.9)
                        {
                            s4 = paste0(loc1,loc1)
                        }else if((n2/(n1+n2)) >= 0.9)
                        {
                            s4 = paste0(loc2,loc2)
                        }else
                        {
                            s4 = paste0(loc1,loc2)
                        }
                        
                    }else if(nchar(s3_unique) == 1)
                    {
                        s4 = paste0(s3_unique,s3_unique)
                    }
                    syn_geno[x,count_col] = s4

                }
                count_col = count_col + 1
            }
            rownames(syn_geno) = rownames(train_geno[[1]])
            colnames(syn_geno) = sort(unique(l1))
            a1 = rrplub_geno(syn_geno)
            ans = kinship.BLUP(y=as.vector(pheno),G.train=a1,G.pred=a2,K.method="GAUSS")
            TrainGeno = syn_geno
            TrainPheno = pheno
        }
        else
        {
            stop("ERROR: incorrect level")
        }
        id = names(ans$g.pred)
        value = ans$g.pred
        predict_pheno = as.data.frame(cbind(id,value))
    }else if(M == "parental_contribution")
    {
        predict_pheno = parental_contribution(pop=predict_geno_pop,parental_genotypes=train_geno,parental_phenotypes=train_pheno)
        TrainGeno = train_geno
        TrainPheno = train_pheno
    }
    return(list(training_sets = list(train_geno = TrainGeno, train_pheno = TrainPheno),phenotypes = predict_pheno))
}

