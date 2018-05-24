create_marker_set <- function(genotypes,N,MAF,no_qtn,print,biased_selection,gen2phy,parental_genotypes)
{
    if(missing(genotypes))
    {
        stop("ERROR: input genotypes should be provided")
    }
    if(missing(MAF) || is.na(MAF))
    {
        maf = 0
    }else
    {
        maf = MAF
    }
    if(missing(N))
    {N = NA}else
    {N = N}
    if(missing(no_qtn))
    {no_qtn = NA}else
    {no_qtn = no_qtn}
    if(missing(print))
    {
        print = TRUE
    }else
    {
        print = print
    }
    if(missing(biased_selection))
    {
        biased_selection = FALSE
    }else
    {
        biased_selection = biased_selection
    }
    if(nrow(genotypes[[1]]) > 100000)
    {
        n1 = sort(sample(1:nrow(genotypes[[1]]),100000))
        A = genotypes[[1]][n1,]
        B = as.data.frame(genotypes[[2]][n1,])
        colnames(B) = "gen_loc"
        rownames(B) = rownames(genotypes[[2]])[n1]
    }else
    {
        A = genotypes[[1]]
        B = genotypes[[2]]
    }
    if(is.character(no_qtn))
    {
        A = A[-(which(rownames(A) %in% no_qtn)),]
        b = which(rownames(B)%in%no_qtn)
        B1 = rownames(B)[-b]
        B2 = B[-b,]
        B = data.frame(B2)
        rownames(B) = B1
    } 
    if(!is.na(N))
    {    
        if(biased_selection == FALSE)
        {
            C = NA
            xx = apply(genotypes[[1]],1,paste0,collapse="")    
            for (i in 1:nrow(genotypes[[1]]))
            {
                xx2 = rawToChar(unique(charToRaw(xx[i])))
                if(nchar(xx2) == 1)
                {
                    xx2 = paste0(xx2,xx2)
                }
	            snp1 = strsplit(xx2,'')[[1]][1]
                snp2 = strsplit(xx2,'')[[1]][2]
                xx_1 = gsub(snp2,'',xx[i])
                xx_2 = gsub(snp1,'',xx[i])
                if(nchar(xx_1)/nchar(xx[i]) >= maf || nchar(xx_2)/nchar(xx[i]) >= maf)
                {
                    C = c(C,names(xx[i]))
                }
            }
            C = C[2:length(C)]
            if(N >= length(C))
            {
                if(print == TRUE){print(paste0("WARNING: all QTNs with maf >= ",MAF," were returned"))}
            }else
            {
                C = C[sort(sample(1:length(C),N))]
                
            }
        }
        else
        {
            C = as.character(random_qtn_assign(qtn = N,gen2phy=gen2phy,biased_selection=TRUE,parental_genotypes=parental_genotypes,min_qtn_freq=maf,dominant="FALSE",effect_distribution="equal")[,1])
        }
        A2 = A[which(rownames(A) %in% C),]
        B2 = as.data.frame(B[which(rownames(B) %in% C),])
        rownames(B2) = rownames(B)[which(rownames(B) %in% C)]
        colnames(B2) = "gen_loc"
    }else
    {
        stop("ERROR: you have to assigned number of markers to be selected")
    }  
    return(list(genotypes=A2,gen2phy=B2))
}

