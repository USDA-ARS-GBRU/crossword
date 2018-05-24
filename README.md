# crossword
A data-driven simulation language for the design of genetic-mapping experiments and breeding strategies

Installation:
1. if simcross is not installed, we suggest to install it in advance
    devtools::install_github("kbroman/simcross")
2. installing crossword          
    devtools::install_github("USDA-ARS-GBRU/crossword")
3. loading crossword
    library('crossword')

Citation: The publication will be published soon

Implementation:

1. The manual is available at Wiki page.

2. for simple example running:
  
  script_file = "crossword_script_input_file.script"
  script_file = paste0(system.file("extdata",package="crossword"),"/",script_file)
  run_pipeline(script_file)

3. for batch running:
  script_file2 = "test2.script"
  script_file2 = paste0(system.file("extdata",package="crossword"),"/",script_file2)
  list2 = "list2"
  list2 = paste0(system.file("extdata",package="crossword"),"/",list2)
  run_batch(script_file2,list2,run=TRUE)
4. for running using GUI:
    crossword_gui()

5. for running using Rscript:
  Rscript ./crossword.R tutoria
  Rscript ./crossword.R script_file

