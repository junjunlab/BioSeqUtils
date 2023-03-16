#' GenomeGTF
#'
#' @author JunZhang
#'
#' @slot gtf GRanges. The GTF/GFF annotation file.
#' @slot genome DNAStringSet. Genome sequence.
#' @slot representTrans data.frame. longest transcript information.
#' @slot intron data.frame. The intron data.
#' @slot gtfPath character. Path for annnotation file.
#' @slot genomePath character. Path for genome file.
#' @slot ranges GroupedIRanges. Null.
#'
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom progress progress_bar
#' @importFrom methods as
#' @return return a GenomeGTF class.
#' @export
setClass(Class = "GenomeGTF",
         slots = c(gtf = "GRanges",
                   gtfPath = "character",
                   genome = "DNAStringSet",
                   genomePath = "character",
                   representTrans = "data.frame",
                   intron = "data.frame",
                   ranges = "GroupedIRanges"),
         prototype = list(gtf = NULL,
                          genome = NULL,
                          gtfPath = NULL,
                          genomePath = NULL),
         contains = c("DNAStringSet"))



#' loadGenomeGTF
#' Load genome annotation information from a GTF or GFF file and a genome fasta file.
#'
#' @author JunZhang
#'
#' @param gtfPath The path to the GTF or GFF annotation file. Default is NULL.
#' @param genomePath The path to the genome fasta file. Default is NULL.
#' @param format The annotation file format ("gtf" or "gff"). Default is "gtf".
#' @param filterProtein Whether to filter protein-coding genes. Default is FALSE.
#'
#' @return A GenomeGTF object containing genome annotation information.
#' @examples
#' \dontrun{
#' # make object
#' mytest <- loadGenomeGTF(gtfPath = "hg38.ncbiRefSeq.gtf.gz",
#'                         genomePath = "hg38.fa.gz")
#' }
#' @export
loadGenomeGTF <- function(gtfPath = NULL,
                          genomePath = NULL,
                          format = "gtf",
                          filterProtein = FALSE){
  # load gtf
  ncbiRefSeq <- rtracklayer::import.gff(gtfPath,format = format)
  # as.data.frame()

  # load genome sequences
  if(!is.null(genomePath)){
    myFASTA <- Biostrings::readDNAStringSet(genomePath,format = "fasta")
  }else{
    genomePath = ""
    myFASTA = Biostrings::DNAStringSet(NULL)
  }

  # whether get protein
  if(filterProtein == TRUE){
    protein <- ncbiRefSeq[which(ncbiRefSeq$type == "CDS"),]
    gtf <- ncbiRefSeq[which(ncbiRefSeq$gene_id %in% unique(protein$gene_id)),]
  }else{
    gtf <- ncbiRefSeq
  }

  # # extarct all genes information from GTF using julia
  # if(runJulia == TRUE){
  #   # export gtf to gff3
  #   rtracklayer::export(gtf, "GenomeGTFobject.gff3", format = "gff3")
  #
  #   # check julia
  #   if(TRUE){
  #     Julia <- JuliaCall::julia_setup(JULIA_HOME = JuliaHome,installJulia = installJulia)
  #
  #     script_path <- paste0('include("',
  #                           system.file("extdata", "filterRepTrans.jl", package = "BioSeqUtils"),
  #                           '")',collapse = "")
  #
  #     # filterRepTrans <- JuliaCall::julia_eval('include("./R/filterRepTrans.jl")')
  #     filterRepTrans <- JuliaCall::julia_eval(script_path)
  #
  #     # excute function
  #     filterRepTrans("GenomeGTFobject.gff3","allTranscriptInfo.txt",sep)
  #
  #     # load results
  #     representTrans <- utils::read.delim("allTranscriptInfo.txt",sep = "\t",header = FALSE)
  #
  #     # add name
  #     colnames(representTrans) <- c("gname","gid","tid","cdsst","cdsed","tlen","tname","cdslen")
  #     representTrans <- representTrans %>%
  #       dplyr::mutate(cdsed = ifelse(cdslen == 0,tlen,cdsed)) %>%
  #       dplyr::mutate(tname = ifelse(cdslen == 0,paste(tname,cdsst,cdsed,tlen,"NC",sep = sep),
  #                                    paste(tname,cdsst,cdsed,tlen,"CD",sep = sep))
  #                     )
  #
  #   }else{
  #     message("Could not find Julia, please make sure Julia has been installed.")
  #     representTrans <- data.frame()
  #   }
  # }else{
  #   representTrans <- data.frame()
  # }


  # initialize
  object <-
  methods::new("GenomeGTF",
               gtf = gtf,
               gtfPath = gtfPath,
               genome =  myFASTA,
               genomePath = genomePath)

  myShow(object)
  return(object)
}
