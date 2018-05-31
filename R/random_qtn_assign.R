random_qtn_assign <- function(qtn,gen2phy,biased_selection,parental_genotypes,min_qtn_freq,dominant,effect_distribution,highest_P,lowest_P,high_to_low_percentage)
{
    if(missing(qtn))
    {
        stop("ERROR: qtn should be provided")
    }
    if(missing(gen2phy))
    {
        print("WARNING: gen2phy is missing, biased_selection can not be applied")
    }
    if(missing(biased_selection))
    {
        biased_selection = FALSE
    }else
    {
        biased_selection = biased_selection
    }
    if(missing(parental_genotypes))
    {
        stop("ERROR: parental_genotypes should be provided")
    }
    if(missing(min_qtn_freq) || is.na(min_qtn_freq))
    {
        min_qtn_freq = 0
    }else
    {
        min_qtn_freq = min_qtn_freq
    }
    if(missing(dominant))
    {
        dominant = FALSE
    }else
    {
        dominant = dominant
    }
    if(missing(effect_distribution))
    {
        effect_distribution = "equal"
    }else if (effect_distribution == "gamma" || effect_distribution == "normal")
    {
        effect_distribution = effect_distribution
    }
    if(missing(highest_P))
    {
        highest_P = NA
    }else
    {
        highest_P=highest_P
    }
    if(missing(lowest_P))
    {
        lowest_P = NA
    }else
    {
        lowest_P=lowest_P
    }
    if((is.na(highest_P) && !is.na(lowest_P)) || !is.na(highest_P) && is.na(lowest_P))
    {
        stop("ERROR, both highest and lowest parentes should be provided together or missed together")
    }
    if(missing(high_to_low_percentage))
    {
        high_to_low_percentage = 0
    }else
    {
        high_to_low_percentage = high_to_low_percentage
    }
    po2 = parental_genotypes$genotypes
    if(nrow(po2) > 10000)
    {
        po2 = po2[sample(1:nrow(po2),10000),]
    }
    ######################################################################################################################
    ## in case of the user needs to apply biased selection toward the regions with high gene density vs non-coding regions
    FREQ=data.frame()
    if(biased_selection == TRUE)
    {
        for (i in 1:nrow(po2))
        {
	        id = rownames(po2)[i]
	        if (!missing(biased_selection) && !missing(gen2phy) && biased_selection == TRUE)
	        {
	            chr = gsub("^(.*)\\_.*.*", "\\1",id)
	            loc = as.numeric(gsub(".*_(.*)\\.*", "\\1",id))
	            freq2 = gen2phy[[chr]]$data[which(gen2phy[[chr]]$data$window_end > loc)[1],4]
	            FREQ[i,1] = i
	            FREQ[i,2] = freq2 
	        }
	    }
	}
	select = NA
    if(!is.na(highest_P) & !is.na(lowest_P))
    {
        h = highest_P
        l = lowest_P
        po = parental_genotypes$genotypes[,c(h,l)]
        po[po[,1]==po[,2],] = NA
    }
    count = 1
    while (TRUE)
    {
        if (!missing(biased_selection) && !missing(gen2phy) && biased_selection == TRUE)
        {
	        x = sample(FREQ$V1,1,prob=FREQ$V2)
        }else 
        {
	        x = sample(1:nrow(po2),1)
        }
        len2 = (length(po2[x,])*2)
        xx = apply(as.data.frame(po2[x,]),2,paste0, collapse="")
        xx2 = rawToChar(unique(charToRaw(xx)))
        snp1 = strsplit(xx2,'')[[1]][1]
        snp2 = strsplit(xx2,'')[[1]][2]
        freq = ((nchar(gsub(snp2,"",xx))) / len2)
        freq2 = ((nchar(gsub(snp1,"",xx))) / len2)
        if (freq > min_qtn_freq && freq2 > min_qtn_freq && !(x %in% select) && length(unique(po2[x,])) > 1 && is.na(highest_P) && is.na(lowest_P))
        {
	        select = c(select,x)
        }else if (!is.na(highest_P) && !is.na(lowest_P) && !is.na(po[x,1]) && !(x %in% select))
        {
            select = c(select,x)
        }
        if (length(select)>(qtn))
        {
	        break
        }
        count = count+1
        if(count == (qtn*1000))
        {
            stop(paste0("ERROR: ",(qtn*1000)," trials were done without capturing requested number of QTNs"))
        }
    }
    select = select[2:length(select)]
	po3 = po2[sort(select),]	
	SNPs = data.frame()
    po4 = matrix(ncol = ncol(po3),nrow=nrow(po3))
    Eff = data.frame()
    for(i in 1:nrow(po3))
    {
	    id = rownames(po3)[i]
	    xx = apply(as.data.frame(po3[i,]),2,paste0, collapse="")
	    xx2 = rawToChar(unique(charToRaw(xx)))
	    snp1 = strsplit(xx2,'')[[1]][1]
        snp2 = strsplit(xx2,'')[[1]][2]
	    SNPs[i,1] = id
	    SNPs[i,2] = snp1
	    SNPs[i,3] = snp2
	    eff = sample(c(1,-1),1)
	    Eff[i,1] = id
	    if(eff==1)
	    {
	        Eff[i,2] = paste0(snp1,snp1)
	    }else
	    {
	        Eff[i,2] = paste0(snp1,snp1)
	    }
	    Eff[i,3] = eff
	    po4[i,which(po3[i,]==paste0(snp1,snp1))] = eff
	    po4[i,which(po3[i,]==paste0(snp2,snp2))] = 0-eff
	    if (missing(dominant) || dominant == FALSE)
	    {
		    po4[i,which(po3[i,]==paste0(snp1,snp2))] = 0
		    po4[i,which(po3[i,]==paste0(snp2,snp1))] = 0

	    }else
	    {
		    po4[i,which(po3[i,]==paste0(snp1,snp2))] = 1
		    po4[i,which(po3[i,]==paste0(snp2,snp1))] = 1
	    }
	}
    rownames(po4) = rownames(po3)
    colnames(po4) = colnames(po3)
    colnames(SNPs) = c('id','A','B')
    ##calculating the effects
    if(effect_distribution == "equal")
    {
        eff <- rep(1,qtn)/qtn
    }else if(effect_distribution == "gamma")
    {
        geff = rgamma(qtn,5,0.2)
        eff = geff/sum(geff)
    }
    rownames(Eff) = Eff[,1]
    colnames(Eff) = c("QTN","ref","effect")
    Eff$effect = eff * Eff$effect
    if (!is.na(highest_P) && !is.na(lowest_P))
    {
        Eff$effect = abs(as.numeric(as.character(Eff$effect)))
        Eff$ref = po2[as.character(Eff$QTN),h]
        if(high_to_low_percentage > 0)
        {
            ind1 = sample(1:nrow(Eff),ceiling(high_to_low_percentage * qtn))
            Eff$effect[ind1]  = Eff$effect[ind1] * -1
        }
    }
    return(Eff)        
}

