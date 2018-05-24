pick_individual <- function(pop,cross,family,individual)
{
    if(missing(pop))
    {
        stop("ERROR: population object should be supplied")
    }
    if(missing(cross))
    {
        cross = NA
    }else
    {
        cross = cross
    }
    if(missing(family))
    {
        family = NA
    }else
    {
        family = family
    }
    if(missing(individual))
    {
        stop("ERROR: an individual should be provided")
    }else
    {
        individual = individual
    }
    a = pop$scheme
    if(!is.na(cross))
    {
        Cross=unique(a$C)
        Cross = as.numeric(as.character(Cross[cross]))
        b = a[a$C == Cross,]
        if(length(unique(as.numeric(as.character(a$C)))) == 1)
        {
            print("WARNING: you selected 1 cross from a population contains 1 cross")
        }
        a=b
    }
    if(!is.na(family))
    {
        Family=unique(a$S)
        Family = as.numeric(as.character(Family[family]))
        b = a[a$S == Family,]
        if(length(unique(as.numeric(as.character(a$S)))) == 1)
        {
            print("WARNING: you selected 1 family from a population contains 1 family")
        }
        a=b
    }
    if(is.na(individual))
    {
        stop("ERROR: individual should be provided")
    }else
    {
        Individual=unique(a$S)
        Individual = as.numeric(as.character(Individual[individual]))
        a = a[a$N==individual,]
    }
    id=as.character(a$id)
    out = select_haplotype(pop,id)
    return(out)
}

