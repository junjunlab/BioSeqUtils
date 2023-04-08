#' Load BED files and create a data frame
#'
#' This function takes one or more BED files and reads them into a data frame
#' containing columns for chromosome, start position, end position, sample name,
#' and sample number.
#'
#' @param bed_file A character vector of one or more file paths to BED files.
#' @param file_name A character vector of file name for re-assign.
#'
#' @return A data frame containing columns for chromosome, start position, end
#' position, sample name, and sample number.
#'
#' @importFrom plyr ldply
#' @importFrom rtracklayer import.bed
#' @importFrom dplyr select
#'
#' @examples
#' \dontrun{bedfiles <- c("sample1.bed", "sample2.bed")
#' bed_df <- loadBed(bedfile = bedfiles)}
#'
#' @export
loadBed <- function(bed_file = NULL,file_name = NULL){
  plyr::ldply(1:length(bed_file),function(x){
    # tmp <- read.table(bedfile[x],sep = "\t")[,1:3]
    tmp <- rtracklayer::import.bed(bed_file[x]) %>% data.frame() %>%
      dplyr::select(seqnames,start,end)

    # colnames(tmp) <- c("seqnames","start","end")
    # add name
    if(is.null(file_name)){
      tmp$sampleName <- strsplit(bed_file[x],split = ".bed") %>% unlist()
    }else{
      tmp$sampleName <- file_name[x]
    }

    # add sn
    tmp$y <- x
    return(tmp)
  }) -> bed_df
}


#' Load BigWig files and filter chromosomes
#'
#' This function reads in one or multiple BigWig files and filters them by
#' chromosome. The output is a data frame that contains the data from the BigWig
#' file(s) and the sample name.
#'
#' @param bw_file character vector of file path to BigWig files to be loaded
#' @param file_name character vector of sample name(s) corresponding to the BigWig
#' file(s). If NULL, the sample name will be extracted from the file path.
#' @param chrom character vector of chromosome name(s) to be included in the
#' output. If NULL, all chromosomes are included.
#'
#' @return A data frame that contains the filtered data from the BigWig file(s)
#' and the sample name.
#'
#' @examples
#' \dontrun{
#' # Load all chromosomes from one BigWig file with a specified sample name
#' loadBigWig(bw_file = "/path/to/file.bw", file_name = "sample1")
#'
#' # Load only specific chromosomes from multiple BigWig files with default sample names
#' loadBigWig(bw_file = c("/path/to/file1.bw", "/path/to/file2.bw"), chrom = c("chr1", "chr2"))
#' }
#'
#' @export
loadBigWig <- function(bw_file = NULL,file_name = NULL,chrom = NULL){
  # loop read bed
  plyr::ldply(1:length(bw_file),function(x){
    tmp <- rtracklayer::import.bw(bw_file[x])

    # filter chromosome
    if(!is.null(chrom)){
      tmp <- data.frame(tmp %>% plyranges::filter(seqnames %in% chrom) %>%
                          data.frame() %>%
                          plyranges::select(-width,-strand))
    }else{
      tmp <- data.frame(tmp)
    }

    # sampe name
    if(is.null(file_name)){
      spt <- strsplit(bw_file[x],split = "/|.bw|.bigwig") %>% unlist()
      sname <- spt[length(spt)]
    }else{
      sname <- file_name[x]
    }

    # add name
    tmp$fileName <- sname

    return(tmp)
  }) -> bWData
  return(bWData)
}
