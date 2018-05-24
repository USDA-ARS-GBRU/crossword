count_chr_length <- function(fa,output)
{
    require("Biostrings")
    s = readDNAStringSet(fa)
    A = data.frame()
    count = 1
    for (i in names(s))
    {
        id = i
        len = length(s[[i]])
        A[count,1] = id
        A[count,2] = len
        count = count + 1
    }
    write.table(A,output,quote = FALSE, sep = "\t",row.names=FALSE,col.names=FALSE)
}

