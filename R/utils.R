globalVariables(c("all_karyotype","ce_end","ce_start","aes", "arrow", "chr",
                  "element_blank", "element_text", "end.x", "end.y", "facet_wrap",
                  "geom_polygon", "geom_rect", "geom_text", "karyotype_df",
                  "label", "numGenes", "scale_fill_gradient", "segy", "start.x",
                  "theme", "theme_bw", "unit", "x", "xlab", "y", "ylab",
                  "PositionIdentity", "StatIdentity", "layer","density"))

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
    seqlen <- Biostrings::width(object@genome[x]) %>% as.character()
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
#' @param gtf A character which contains the path of a gtf file.
#' @param feature The feature to be extracted for density calculation, default "transcript".
#' @param step A numerical value which represents the length of each segmentation,
#' default is 5,000,000bp.
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
                          gtf = NULL,
                          feature = "transcript",
                          step = 5*10^6){
  # ============================================
  # 1_get intervals
  plyr::ldply(1:nrow(karyotype),function(x){
    tmp <- karyotype[x,]
    res <- calcuIntervals(tmp$chr,tmp$end,step) %>%
      mutate(type = tmp$type)
  }) -> all_intervals

  # ============================================
  # 2_get gene densitys
  if(is.character(gtf)){
    gene <- data.frame(rtracklayer::import.gff(gtf,format = "gtf"))
  }else{
    gene <- gtf
  }

  gene <- gene[which(gene$type == feature & gene$seqnames %in% unique(karyotype$chr)),] %>%
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
  gene_density$density <- gene_density$numGenes/sum(gene_density$numGenes)
  return(gene_density)
}


#' Convert and transform coordinates of genomic features in a GTF file
#'
#' This function takes a GTF file as input and converts the coordinates of
#' genomic features to transcript coordinates. It also reassigns the start and
#' end coordinates for exons, as well as for coding regions (CDS).
#'
#' @param gtf_file A data frame object with at least the columns "transcript_id",
#' "type", "strand", "start", "end", and "width".
#'
#' @return A data frame with the same columns as \code{gtf_file}, where the
#' coordinates of features are in transcript space.
#'
#' @examples
#' \dontrun{
#' # Run function
#' transCoordTransform(gtf_file = example_gtf)
#' }
#'
#' @export
transCoordTransform <- function(gtf_file = NULL){
  tid <- unique(gtf_file$transcript_id)

  # x = 1
  plyr::ldply(seq_along(tid),
              .parallel = TRUE,
              function(x){
                tmp <- gtf_file[which(gtf_file$transcript_id %in% tid[x]),]

                # params
                feature <- unique(tmp$type)
                strand <- unique(tmp$strand)
                trans_len <- sum(tmp[which(tmp$type == "exon"),"width"])

                # order
                if(strand == "+"){
                  tmp <- tmp %>% arrange(start,end)
                }else{
                  tmp <- tmp %>% arrange(-start,-end)
                }

                # re-assign transcript coord
                transcript_tmp <- tmp[which(tmp$type == "transcript"),] %>%
                  mutate(start = 1,end = trans_len)

                # re-assign coord for "exon" and "CDS/3UTR/5UTR"
                # for exon
                exon_tmp <- tmp[which(tmp$type == "exon"),]
                exon_len <- cumsum(exon_tmp$width)
                new_start <- c(0,exon_len[1:(length(exon_len) - 1)]) + 1

                # check whether only one feature
                if(nrow(exon_tmp) == 1){
                  exon_tmp$start <- new_start[1]
                  exon_tmp$end <- exon_len
                }else{
                  exon_tmp$start <- new_start
                  exon_tmp$end <- exon_len
                }

                # for CDS
                if("CDS" %in% feature){
                  coding_tmp <- tmp[which(tmp$type %in% c("5UTR","five_prime_utr",
                                                          "CDS",
                                                          "3UTR","three_prime_utr")),]
                  exon_len <- cumsum(coding_tmp$width)
                  new_start <- c(0,exon_len[1:(length(exon_len) - 1)]) + 1

                  if(nrow(exon_tmp) == 1){
                    coding_tmp$start <- new_start[1]
                    coding_tmp$end <- exon_len
                  }else{
                    coding_tmp$start <- new_start
                    coding_tmp$end <- exon_len
                  }
                }else{
                  coding_tmp <- NULL
                }

                # combine data
                fina_res <- rbind(transcript_tmp,exon_tmp,coding_tmp)

                return(fina_res)
              }) -> trans_gtf
  return(trans_gtf)
}



# =========================================================================================
# test data
# =========================================================================================

#' This is a test data for this package
#' test data describtion
#'
#' @name all_karyotype
#' @docType data
#' @author JunZhang
"all_karyotype"
