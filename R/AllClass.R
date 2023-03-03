#' GenomeGTF
#'
#' @author JunZhang
#' @slot gtf GRanges. The GTF/GFF annotation file.
#' @slot genome DNAStringSet. Genome sequence.
#' @slot representTrans data.frame. longest transcript information.
#' @slot intron data.frame. The intron data.
#' @slot ranges GroupedIRanges. Null.
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom progress progress_bar
#' @importFrom methods as
#' @return return a GenomeGTF class.
#' @export
setClass(Class = "GenomeGTF",
         slots = c(gtf = "GRanges",
                   genome = "DNAStringSet",
                   representTrans = "data.frame",
                   intron = "data.frame",
                   ranges = "GroupedIRanges"),
         prototype = list(gtf = NULL,
                          genome = NULL),
         contains = c("DNAStringSet"))



#' loadGenomeGTF
#' Load genome annotation information from a GTF or GFF file and a genome fasta file.
#'
#' @author JunZhang
#' @param gtfPath The path to the GTF or GFF annotation file.
#' @param genomePath The path to the genome fasta file.
#' @param format The annotation file format ("gtf" or "gff"), default is "gtf".
#' @param filterProtein Whether to filter protein-coding genes, default is FALSE.
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
  myFASTA <- Biostrings::readDNAStringSet(genomePath,format = "fasta")

  # whether get protein
  if(filterProtein == TRUE){
    protein <- ncbiRefSeq[which(ncbiRefSeq$type == "CDS"),]
    gtf <- ncbiRefSeq[which(ncbiRefSeq$gene_id %in% unique(protein$gene_id)),]
  }else{
    gtf <- ncbiRefSeq
  }

  # filter all introns
  # allintron <- getIntronInfo(mytest,transId = unique(gtf$transcript_id))

  # initialize
  object <-
  methods::new("GenomeGTF",
               gtf = gtf,
               genome =  myFASTA)

  myShow(object)
  return(object)
}
