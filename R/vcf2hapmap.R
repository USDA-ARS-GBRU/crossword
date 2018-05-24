vcf2hapmap <- function(input,output)
{
    if(missing(input))
    {
        stop("ERROR: input should be provided")
    }
    if(missing(output))
    {
        output = paste0(input,".hapmap")
    }
    asd = get_parental_genotypes(input)
    genotypes_out(asd,output)
}

