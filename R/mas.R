mas <- function(pop,marker,level,parental_genotypes)
{
    if(missing(pop))
    {
        stop("ERROR: a hapltype object should be provided")
    }
    if(missing(marker))
    {
        stop("ERROR: a marker object should be provided")
    }
    if(missing(level))
    {
        stop("ERROR: the level should be assigned")
    }
     if(missing(parental_genotypes))
    {
        stop("ERROR: parental_genotypes should be provided")
    }
    scheme = pop[[1]]
    Pop = pop
    scheme = scheme[scheme$gen==scheme$gen[nrow(scheme)],]
    Pop = select_haplotype(Pop,scheme$id)
    geno = get_genotypes(parental_genotypes=parental_genotypes,pop=Pop)[[1]]
    geno = geno[,scheme$id]
    IDS = NA
    if (level=="individual" || level == "cross" || level == "family" || level == "population")
    { 
        l1 = get_level(pop$scheme,level)
        for (i in 1:length(marker))
        {
            m = marker[[i]][2]
            maf = as.numeric(marker[[i]][3])
            for (i2 in sort(unique(l1)))
            {
                ids = as.character(scheme[l1 == i2,1])
                g2 = geno[marker[[i]][1],ids]
                n1 = nchar(paste(g2,collapse=""))
                n2 = nchar(gsub(m,"",paste(g2,collapse="")))
                n3 = (n1-n2)/n1
                if(n3 >= maf)
                {
                    IDS = c(IDS,ids) 
                }
            }
        }
        IDS = sort(unique(IDS[!is.na(IDS)]))
        Pop = select_haplotype(Pop,IDS)
        if(length(IDS) == 0)
        {
            print("WARNING: empty output, no individual has the criteria")
        }
    }else
    {
        stop("ERROR: you entered an incorrect level")
    }
    return(Pop)
}

