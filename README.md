# crossword
A data-driven simulation language for the design of genetic-mapping experiments and breeding strategies

## Installation:
0. Open a command prompt and run R by typing "R[enter]".  Unless otherwise indicated, the following commands should be typed at the R prompt. "devtools" is required to install any package from github. If you do not have it, it can be installed by the following command:

        install.packages("devtools")

1. If simcross is not installed, we suggest to install it in advance.  
    
        devtools::install_github("kbroman/simcross")

2. Install crossword.         
    
        devtools::install_github("USDA-ARS-GBRU/crossword")

## Running crossword:

The full manual is available at Wiki page in the tab above.

### For running using GUI:
            
            library('crossword')
            crossword_gui()
           
The x-window environment should open and commands can be entered accordingly.

### For running directly from crossword script:
  
          library('crossword')
          script_file = paste0(system.file("extdata",package="crossword"),"/crossword_script_input_file.script")
          run_pipeline(script_file)

### For running multiple simulations using a range of parameters:

See the listOfParameters.txt example file.  The variable names in listOfParameters.txt should match those use in the crossword script.  We recommend using all-caps for these names in order to be clear.

          script_file2 = paste0(system.file("extdata",package="crossword"),"/test2.script")
          list = paste0(system.file("extdata",package="crossword"),"/listOfParameters.txt")
          run_batch(script_file2,list2,run=TRUE)

### For running using Rscript (outside R):

These commands are run at the command prompt (not in R, as the above).  Copy "crossword.R" and the script file from the installed library to the local location and modify as needed.
          
          Rscript ./crossword.R tutorial
          Rscript ./crossword.R crossword_script_input_file.script

Auxiliary functions can be run through Rscript by passing the function name then the function's arguments in their order.
            
            Rscript ./crossword.R vcf2hapmap peanut.vcf peanut.hapmap

## Citation: 
The publication has been submitted and will be linked here when accepted.
