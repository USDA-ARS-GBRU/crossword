args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0)
{
    stop("ERROR: an argument required")
}else
{
    library(crossword)
}
if(args[1] == "tutorial")
{
    script_file = "crossword_script_input_file.script"
    script_file = paste0(system.file("extdata",package="crossword"),"/",script_file)
    ptm <- proc.time()
    run_pipeline(script_file)
    (proc.time() - ptm)[1]/60
}else if(args[1] == "filter_parents")
{
    filter_parents(args[2],args[3],args[4])
}else if(args[1] == "count_chr_length")
{
    count_chr_length(args[2],args[3])
}else if(args[1] == "vcf2hapmap")
{
    vcf2hapmap(args[2],args[3])
}else if(args[1] == "create_genome")
{
    create_genome(args[2],args[3],args[4],eval(parse(text=args[5])),as.numeric(args[6]))
}else if(args[1] == "simulate_fastq")
{
    simulate_fastq(as.logical(args[2]),args[3],args[4],as.numeric(args[5]),as.numeric(args[6]),args[7])
}else
{
    script_file = args[1]
    ptm <- proc.time()
    run_pipeline(script_file)
    (proc.time() - ptm)[1]/60
}
