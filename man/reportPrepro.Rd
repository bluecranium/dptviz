\name{reportPrepro}
\alias{reportPrepro}
\title{Summary plots of data preprocessed in bins by iDPT}
\usage{
  reportPrepro(dptpath,
    chrs = c(sapply(seq(1, 22), function(x) paste("chr", x, sep = "")), "chrX", "chrY"),
    bsgenome, file = NULL, ncore = 1)
}
\arguments{
  \item{dptpath}{Path to a folder containing an iDPT run
  (the folder than contains the "DPT_ws" folder)}

  \item{chrs}{A vector of chromosome names to use (by
  default, the autosomes and sex chromosomes from hg19)}

  \item{bsgenome}{A BSgenome object for the reference
  genome used in the alignments contained in the input BAM
  files given to iDPT (i.e. for UCSC's hg19, you would load
  "library(BSgenome.Hsapiens.UCSC.hg19)" and give
  BSgenome.Hsapiens.UCSC.hg19 to this argument). This is
  needed to compute CpG density in the windows used by
  iDPT.}

  \item{file}{Filename to save the PDF report to (if NULL,
  will use the name of the folder as a prefix and save
  "prefix.reportPrepro.pdf" in the current working
  directory")}

  \item{ncore}{Number of concurrent threads to use}
}
\description{
  iDPT performs a preprocessing stage to convert a set of
  alignments in a BAM file to average coverage values in
  windows of a user-specified size genome-wide. This
  function saves a PDF file containing plots that summarize
  the distribution of these average coverage values and
  their relationship to CpG density on a specified
  reference genome.
}

