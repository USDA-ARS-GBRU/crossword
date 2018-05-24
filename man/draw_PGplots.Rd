% File library/man/draw_PGplots.Rd
\name{draw_PGplots}
\alias{draw_PGplots}
\title{draw_PGplots}
\description{
draw_PGplots function draws the qplots of physical to genomical locations for all chrosomes in jpeg format into a folder created by the function
}
\usage{
draw_PGplots(input,output_folder,im_type)
}
\arguments{
\item{input}{Required: the output of physical2genomic function.}
\item{output_folder}{Required: a folder to be created which will be used to store the images.}
\item{im_type}{Optional: the type of images, default is svg}}

\examples{
a1 = physical2genomic("peanut.gff","chr_siz.txt","chr_len.txt",200000)
draw_PGplots(a1,'output_folder_1')}
