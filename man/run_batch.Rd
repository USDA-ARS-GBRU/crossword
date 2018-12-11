% File library/man/run_batch.Rd
\name{run_batch}
\alias{run_batch}
\title{run_batch}
\description{
The function runs a batch crossword scripts.}
\usage{
run_batch(script_file,argument_file,run)
}
\arguments{
\item{script_file}{Required: an input script file}
\item{argument_file}{Required: a file contains the argument description}
\item{run}{Optional: TRUE if the user needs to run the all scripts, default is FALSE}}
\examples{
run_batch("test2.script","listOfParameters",run=TRUE)}
