globalVariables(c("all_karyotype","ce_end","ce_start"))

#' Get chromosome sizes from a genome object
#'
#' This function takes a Genome object and returns a data.frame containing the
#' sizes of the chromosomes.
#'
#' @param object A Genome object
#'
#' @return A data.frame with columns for chromosome name, start and end positions.
#'
#' @importFrom plyr ldply
#' @export
getChromSize <- function(object){
  chr.name <- names(object@genome)
  plyr::ldply(chr.name,function(x){
    seqlen <- nchar(object@genome[x])
    res <- data.frame(chr = x,start = 0,end = seqlen)
    return(res)
  }) -> cl
  return(cl)
}


#' prepareKaryotype
#'
#' This function creates a karyotype data frame for a given chromosome and
#' centromere size.
#'
#' @param chromSize A data frame with columns "chr", "start", and "end" containing
#' the chromosome name, start position, and end position, respectively.
#' @param centromere A data frame with columns "chr", "ce_start", and "ce_end"
#' containing the chromosome name, centromere start position, and centromere end
#' position, respectively. If not supplied, the centromere start and end positions
#' will be set to 0.
#'
#' @return A data frame with columns "chr", "start", "end", "ce_start", and
#' "ce_end" containing the chromosome name, start position, end position,
#' centromere start position, and centromere end position, respectively.
#' @export
prepareKaryotype <- function(chromSize = NULL,
                             centromere = NULL){
  colnames(chromSize) <- c("chr","start","end")
  # check centromere whether supplied
  if(!is.null(centromere)){
    colnames(centromere) <- c("chr","ce_start","ce_end")
    mer <- merge(chromSize,centromere,by = "chr")
  }else{
    mer <- chromSize %>%
      mutate(ce_start = 0,ce_end = 0)
  }
  return(mer)
}


#' Calculate the segmented interval for chromosome
#'
#' Generate segmented intervals for chromosome based on the total length and
#' specified segmentation size.
#'
#' @param chr A character which represents the chromosome of segmented intervals.
#' @param length A numerical value which represents the length of the chromosome.
#' @param step A numerical value which represents the length of each segmentation,
#' default is 1,000,000bp.
#'
#' @return a dataframe which includes three columns representing the chromosome
#' number, starting coordinate, and ending coordinate of each segmented interval.
calcuIntervals <- function(chr,length,step = 2*10^6){
  num.iter <- ceiling(length/step) - 1
  plyr::ldply(0:num.iter, function(x){
    start <- x*step + 1
    if(x == num.iter){
      end <- length
    }else{
      end <- (x+1)*step
    }
    data.frame(chr,start,end)
  }) -> res
  return(res)
}


#' Calculate the centromere data for karyotype
#'
#' Generate data for centromere visualization in karyotype plot.
#'
#' @param karyotype A data frame with four columns：chr, size, ce_start, ce_end.
#' The columns represent the chromosome number, size of the chromosome, start
#' position of the centromere, and end position of the centromere, respectively.
#'
#' @return a data frame which contains three columns: chr, x, and y. x and y
#' represent the coordinate for visualization. The data frame is used for centromere
#' visualization in karyotype plot.
#'
#' @importFrom plyr ldply
#'
#' @export
calcuCentromereData <- function(karyotype = NULL){
  plyr::ldply(1:nrow(karyotype),function(x){
    tmp <- karyotype[x,] %>%
      mutate(ce_mid = (ce_start + ce_end)/2)
    res <- data.frame(chr = tmp$chr,
                      x = c(tmp$ce_start,tmp$ce_start,tmp$ce_mid,tmp$ce_end,tmp$ce_end,
                            tmp$ce_mid,tmp$ce_start),
                      y = c(0,1,0.5,1,0,0.5,0))
    return(res)
  }) -> centromere.triangle
  return(centromere.triangle)
}


#' Get gene density data for karyotype object
#'
#' Get gene density data for a list object of karyotype or a gtf file.
#'
#' @param karyotype A data frame with three columns: chr, start, and end. The
#' columns represent the chromosome number, start position and end position of
#' the chromosome, respectively.
#' @param object A character which contains the path of a gtf file.
#' @param step A numerical value which represents the length of each segmentation,
#' default is 1,000,000bp.
#'
#' @return a data frame which contains four columns: chr, start, end, and numGenes.
#' numGenes represents the number of genes in the corresponding segment.
#'
#' @importFrom plyr ldply
#' @importFrom rtracklayer import.gff
#' @import dplyr
#'
#' @export
getGeneDensiy <- function(karyotype = NULL,
                          object = NULL,
                          step = 2*10^6){
  # ============================================
  # 1_get intervals
  plyr::ldply(1:nrow(all_karyotype),function(x){
    tmp <- all_karyotype[x,]
    res <- calcuIntervals(tmp$chr,tmp$end) %>%
      mutate(type = tmp$type)
  }) -> all_intervals

  # ============================================
  # 2_get gene densitys
  if(is.character(object)){
    gene <- data.frame(rtracklayer::import.gff(object,format = "gtf"))
  }else{
    gene <- data.frame(object@gtf)
  }

  gene <- gene[which(gene$type == "transcript" & gene$seqnames %in% unique(karyotype$chr)),] %>%
    select(seqnames,start,end) %>%
    mutate(mid = (start+end)/2)

  # x = 1
  plyr::ldply(1:nrow(all_intervals),function(x){
    tmp <- all_intervals[x,]
    gene.tmp <- gene[which(gene$seqnames %in% tmp$chr),]
    ng <- table(gene.tmp$mid >= tmp$start & gene.tmp$mid <= tmp$end)[2]
    res <- tmp %>%
      mutate(numGenes = ng)
    res[is.na(res)] <- 0
    return(res)
  }) -> gene_density
}


#' This is a test data for this package
#' test data describtion
#'
#' @name all_karyotype
#' @docType data
#' @author JunZhang
"all_karyotype"
