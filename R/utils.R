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

#' Create Trapezoid Shape
#'
#' This function creates a trapezoid shape with given coordinates.
#'
#' @param xPos A numeric vector of length 4 representing x-coordinates of the
#' trapezoid in order.
#' @param yPos A numeric vector of length 4 representing y-coordinates of the
#' trapezoid in order.
#' @param nDivisions An integer value indicating the number of divisions to
#' create the trapezoid shape. Default is 2000.
#' @param reverse A boolean value indicating whether to create the trapezoid in
#' reverse order or not. Default is FALSE.
#'
#' @return A data frame with id, x and y coordinates.
#'
#' @importFrom plyr ldply
#' @import tidyr
#' @export
createTrapezoid <- function(xPos = NULL,
                            yPos = NULL,
                            nDivisions = 2000,
                            reverse = FALSE){
  n = nDivisions
  # j = 1
  plyr::ldply(1:n,function(j){
    # whether reverse Trapezoid graph
    if(reverse == TRUE){
      x1 = xPos[2];x2 = xPos[1];x3 = xPos[4];x4 = xPos[3]
    }else{
      x1 = xPos[1];x2 = xPos[2];x3 = xPos[3];x4 = xPos[4]
    }
    y1 = yPos[1];y2 = yPos[2];y3 = yPos[3];y4 = yPos[4]

    # calculate positions
    xp <- c(x1 + (x2-x1)/n*(j-1),x1 + (x2-x1)/n*j,x4 - (x4-x3)/n*j,x4 - (x4-x3)/n*(j-1))
    yp <- c(y1 + (y2-y1)/n*(j-1),y1 + (y2-y1)/n*j,y1 + (y2-y1)/n*j,y4 + (y3-y4)/n*(j-1))

    # add id
    res <- data.frame(id = j,x = xp,y = yp)
    return(res)
  }) -> trapezoid.df
  return(trapezoid.df)
}


#' Create segment coordinates
#'
#' This function takes starting and ending coordinates and creates segments in
#' equal divisions or given relative length. It returns a data frame with the x
#' and y coordinates of each segment.
#'
#' @param xPos A vector of length 2 containing the x-coordinates of the starting
#' and ending positions for the segment.
#' @param yPos A vector of length 2 containing the y-coordinates of the starting
#' and ending positions for the segment.
#' @param n_division An integer value indicating the number of segments the line
#' should be divided into.
#' @param rel_len A decimal value indicating the relative length of each segment.
#' @param abs_len A absolute length for generating segments.
#'
#' @return A data frame with the x and y-coordinates of each segment.
#'
#' @examples
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), n_division = 4)
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), rel_len = 0.2)
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), abs_len = 0.1)
#'
#' @export
createSegment <- function(xPos = NULL,
                          yPos = NULL,
                          n_division = NULL,
                          rel_len = NULL,
                          abs_len = NULL){
  if(xPos[1] > xPos[2]){
    x = xPos[2];xend = xPos[1]
  }else{
    x = xPos[1];xend = xPos[2]
  }

  if(yPos[1] > yPos[2]){
    y = yPos[2];yend = yPos[1]
  }else{
    y = yPos[1];yend = yPos[2]
  }

  # get seqs
  if(y == yend){
    seg_len = abs(x - xend)
    if(!is.null(n_division)){
      seqs <- seq(from = x,to = xend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% seqs)) seqs <- c(seqs,xend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = seqs[1:(length(seqs) - 1)],
                      xend = seqs[2:length(seqs)],
                      y = y,
                      yend = y)
  }else if(x == xend){
    seg_len = abs(y - yend)
    if(!is.null(n_division)){
      seqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = y,to = yend,by = by)
      if(!(yend %in% seqs)) seqs <- c(seqs,yend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = x,
                      xend = x,
                      y = seqs[1:(length(seqs) - 1)],
                      yend = seqs[2:length(seqs)])
  }else{
    xseg_len = abs(xend - x)
    yseg_len = abs(yend - y)
    if(!is.null(n_division)){
      xseqs <- seq(from = x,to = xend,length.out = n_division)
      yseqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*xseg_len
      }else{
        by = abs_len
      }
      xseqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% xseqs)) xseqs <- c(xseqs,xend)
      yseqs <- seq(from = y,to = yend,length.out = length(xseqs))
    }
    res <- data.frame(id = 1:(length(xseqs) - 1),
                      x = xseqs[1:(length(xseqs) - 1)],
                      xend = xseqs[2:length(xseqs)],
                      y = yseqs[1:(length(yseqs) - 1)],
                      yend = yseqs[2:length(yseqs)])
  }
  return(res)
}


#' Create a Pair of Polygons
#'
#' This function creates a pair of polygons, each consisting of two line segments,
#' with a shared mid-point and a user-defined start and end point.
#'
#' @param start A numeric value indicating the start point of the polygon.
#' Defaults to NULL.
#' @param end A numeric value indicating the end point of the polygon. Defaults
#' to NULL.
#' @param open A numeric value indicating the openness of the polygon. Defaults
#' to 1.
#' @param npoints An integer value indicating the number of points to be used to
#' construct the polygon. Defaults to 1000.
#' @param reverse A logical value indicating whether to reverse the order of the
#' polygons.
#'
#' @return A data.frame object containing the x and y coordinates of the two
#' polygons and the type of polygon.
#'
#' @examples
#' createPairPolygon(start = 0, end = 10, open = 1, npoints = 1000, reverse = FALSE)
#'
#' @export
createPairPolygon <- function(start = NULL,
                              end = NULL,
                              open = 1,
                              npoints = 1000,
                              reverse = FALSE){
  mid = (start + end)/2
  interval_1 = seq(start,mid,length.out = npoints)

  # check start+end = 0
  if(mid == 0){
    y_interval = abs(seq(0,sum(abs(start) + abs(end))/2,
                         length.out = npoints))
  }else{
    y_interval = seq(0,mid,length.out = npoints)
  }

  ratio = sqrt(abs(max(y_interval))/2)
  x = c(interval_1,rev(interval_1))
  y = c(-(sqrt(abs(rev(y_interval)/2))*open),
        sqrt(abs(y_interval)/2)*open)/ratio


  interval_2 = seq(mid,end,length.out = npoints)
  x1 = c(rev(interval_2),interval_2)
  y1 = y

  if(reverse == TRUE){
    ypos = c(x1,x);xpos = c(y,y1)
  }else{
    xpos = c(x1,x);ypos = c(y,y1)
  }

  res <- data.frame(xpos = xpos,ypos = ypos,
                    type = rep(c(1,2),each = length(x)))
}



#' Add a custom annotation layer to a ggplot
#'
#' This function adds a custom annotation layer to a ggplot using the provided grob object.
#' The annotation can be placed within specified x and y limits.
#'
#' @param grob A grob object that will be used as the custom annotation.
#' @param xmin The minimum x value where the custom annotation will be placed (default is -Inf).
#' @param xmax The maximum x value where the custom annotation will be placed (default is Inf).
#' @param ymin The minimum y value where the custom annotation will be placed (default is -Inf).
#' @param ymax The maximum y value where the custom annotation will be placed (default is Inf).
#' @param data The data that the custom annotation layer will be applied to.
#' @return A ggplot layer with the custom annotation.
#' @export
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob,
                                          xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax))
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
