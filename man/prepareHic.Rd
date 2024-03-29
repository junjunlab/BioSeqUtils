% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importFileFunctions.R
\name{prepareHic}
\alias{prepareHic}
\title{Prepare Hi-C data for analysis}
\usage{
prepareHic(
  hic_path = NULL,
  readHic_params = list(),
  data = NULL,
  file_name = NULL,
  assembly = NULL,
  chrom = NULL,
  chromstart = NULL,
  chromend = NULL,
  resolution = NULL
)
}
\arguments{
\item{hic_path}{character vector specifying the path(s) to Hi-C data file(s).}

\item{readHic_params}{list of parameters to be passed to the
\code{\link[plotgardener]{readHic}} function, which is used to read in Hi-C
data files. See the documentation for \code{\link[plotgardener]{readHic}} for
more information.}

\item{data}{data frame containing Hi-C data. This parameter is used if Hi-C
data has already been loaded and does not need to be read in from file(s).}

\item{file_name}{character vector specifying the file name(s) associated with
the Hi-C data file(s).}

\item{assembly}{character vector specifying the assembly(ies) of the Hi-C data.}

\item{chrom}{character vector specifying the chromosome(s) of the Hi-C data.}

\item{chromstart}{numeric vector specifying the start position(s) of the Hi-C data.}

\item{chromend}{numeric vector specifying the end position(s) of the Hi-C data.}

\item{resolution}{numeric vector specifying the resolution(s) of the Hi-C data.}
}
\value{
a data frame containing the prepared Hi-C data.
}
\description{
This function prepares Hi-C data for downstream analysis by either reading
in Hi-C data files or using pre-loaded data.
}
\examples{
\dontrun{
prepareHic(hic_path = "my_hic_file.hic", assembly = "hg38", chrom = "chr1", resolution = 5000)
}

}
