globalVariables(c("fileName","value","xp","yp",":=",".N"))

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
#' @param format the signal data format, "bw"(default), "wig" or "bedGraph".
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
#' @importFrom fastplyr f_filter f_select
#'
#' @export
loadBigWig <- function(bw_file = NULL,file_name = NULL,chrom = NULL,format = c("bw","bedGraph")){
  format <- match.arg(format,choices = c("bw","bedGraph"))
  # loop read bed
  plyr::ldply(1:length(bw_file),function(x){
    tmp <- rtracklayer::import(bw_file[x],format = format)

    # filter chromosome
    if(!is.null(chrom)){
      # tmp <- data.frame(tmp %>% plyranges::filter(seqnames %in% chrom) %>%
      #                     data.frame() %>%
      #                     plyranges::select(-width,-strand))
      tmp <- data.frame(tmp) %>%
        fastplyr::f_filter(seqnames %in% chrom) %>%
        fastplyr::f_select(-width,-strand)
    }else{
      tmp <- data.frame(tmp) %>% plyranges::select(-width,-strand)
    }

    # sampe name
    if(is.null(file_name)){
      if(format == "bw"){
        fixchar <- "/|.bw|.bigwig"
      }else{
        fixchar <- "/|.bg|.bedgraph"
      }

      spt <- strsplit(bw_file[x],split = fixchar) %>% unlist()
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


#' Prepare Hi-C data for analysis
#'
#' This function prepares Hi-C data for downstream analysis by either reading
#' in Hi-C data files or using pre-loaded data.
#'
#' @param hic_path character vector specifying the path(s) to Hi-C data file(s).
#' @param readHic_params list of parameters to be passed to the
#' \code{\link[plotgardener]{readHic}} function, which is used to read in Hi-C
#' data files. See the documentation for \code{\link[plotgardener]{readHic}} for
#' more information.
#' @param data data frame containing Hi-C data. This parameter is used if Hi-C
#' data has already been loaded and does not need to be read in from file(s).
#' @param file_name character vector specifying the file name(s) associated with
#' the Hi-C data file(s).
#' @param assembly character vector specifying the assembly(ies) of the Hi-C data.
#' @param chrom character vector specifying the chromosome(s) of the Hi-C data.
#' @param chromstart numeric vector specifying the start position(s) of the Hi-C data.
#' @param chromend numeric vector specifying the end position(s) of the Hi-C data.
#' @param resolution numeric vector specifying the resolution(s) of the Hi-C data.
#'
#' @return a data frame containing the prepared Hi-C data.
#'
#' @examples
#' \dontrun{
#' prepareHic(hic_path = "my_hic_file.hic", assembly = "hg38", chrom = "chr1", resolution = 5000)
#' }
#'
#' @importFrom plotgardener readHic
#' @importFrom plyr ldply
#' @import dplyr
#' @import tidyr
#' @export
prepareHic <- function(hic_path = NULL,readHic_params = list(),
                       data = NULL,file_name = NULL,assembly = NULL,
                       chrom = NULL,chromstart = NULL,chromend = NULL,
                       resolution = NULL){
  time_n = max(length(file_name),length(data),length(hic_path))
  if(length(chrom) == 1) chrom = rep(chrom,time_n)
  if(length(resolution) == 1) resolution = rep(resolution,time_n)
  if(length(assembly) == 1) assembly = rep(assembly,time_n)

  if(!is.null(chromstart) & !is.null(chromend)){
    if(length(chromstart) == 1) chromstart = rep(chromstart,time_n)
    if(length(chromend) == 1) chromend = rep(chromend,time_n)
  }

  # read hic data
  if(!is.null(hic_path)){
    if(endsWith(hic_path[1],".hic")){
      # ==============================================================
      # read hic format data
      # x = 1
      hic_data <- lapply(seq_along(hic_path), function(x){
        if(!is.null(chromstart) & !is.null(chromend)){
          readHic_list <- list(file = hic_path[x],
                               chrom = chrom[x],
                               chromstart = chromstart[x],
                               chromend = chromend [x],
                               assembly = assembly[x],
                               resolution = resolution[x])
        }else{
          readHic_list <- list(file = hic_path[x],
                               chrom = chrom[x],
                               assembly = assembly[x],
                               resolution = resolution[x])
        }

        tmp <- do.call(plotgardener::readHic,modifyList(readHic_list,readHic_params))
        return(tmp)
      })
    }else if(endsWith(hic_path[1],".cool")){
      # ==============================================================
      # read cool format data
      if(!is.null(chromstart) & !is.null(chromend)){
        hic_data <- lapply(seq_along(hic_path), function(x){
          # load .cool data
          cool_bedpe <- HiCcompare::cooler2bedpe(path = hic_path[x])$cis[[chrom[x]]]
          # to up-triangle matrix
          sparse <- HiCcompare::cooler2sparse(cool_bedpe)[[chrom[x]]]
          sparse <- sparse[which(sparse$region1 >= chromstart[x] & sparse$region2 >= chromend[x]),]
          # sparse_df <- cbind(seqnames = rep(chrom[x],nrow(sparse)),sparse)
          rhdf5::h5closeAll()
          return(sparse)
        })
      }else{
        hic_data <- lapply(seq_along(hic_path), function(x){
          # load .cool data
          cool_bedpe <- HiCcompare::cooler2bedpe(path = hic_path[x])$cis[[chrom[x]]]
          # to up-triangle matrix
          sparse <- HiCcompare::cooler2sparse(cool_bedpe)[[chrom[x]]]
          # sparse_df <- cbind(seqnames = rep(chrom[x],nrow(sparse)),sparse)
          rhdf5::h5closeAll()
          return(sparse)
        })
      }
    }else{
      message("please supply hic data with .hic or .cool format!")
    }
  }else{
    hic_data <- data
  }

  # x = 1
  plyr::ldply(1:length(hic_data),function(x){
    colnames(hic_data[[x]]) <- c('x','y','count')
    tmp <- getRotatedPolygon(data = hic_data[[x]],rx = 'x',ry = 'y',
                             value = 'count',workers = 1,
                             window = resolution[x])$polygon_coods

    # filename
    if(is.null(file_name)){
      file_name_tmp <- unlist(strsplit(hic_path[x],split = "/|.hic|.cool"))
      file_name_ex <- file_name_tmp[length(file_name_tmp)]
    }else{
      file_name_ex <- file_name[x]
    }

    # formatter data
    tmp$seqnames <- ifelse(endsWith(hic_path[x],".hic"),
                           paste("chr",chrom[x],sep = ""),chrom[x])
    tmp$fileName <- file_name_ex
    tmp <- tmp %>% select(seqnames,xp,yp,value,fileName,id)
    colnames(tmp)[2:5] <- c("start","end","score","fileName")

    return(tmp)
  }) -> hicdf
}


#' Load loop data from bed or bedpe format files
#'
#' This function loads loop data from bed or bedpe format files and returns them
#' as a data frame.
#'
#' @param loop_file a character vector of file names to be loaded.
#' @param file_name a character vector of corresponding names to the file names
#' in \code{loop_file},or \code{NULL} to use the file names themselves as the
#' names of the output data frames.
#' @param sep the separator used in the input files.
#'
#' @return a data frame containing the loop data. The columns of the data frame
#' are as follows:
#' \itemize{
#' \item \code{seqnames}: chromosome names of the interacting regions;
#' \item \code{start}: mid-point of the first region;
#' \item \code{end}: mid-point of the second region;
#' \item \code{score}: a score calculated as the distance between the start and end positions
#' of the interacting regions divided by 100,000;
#' \item \code{fileName}: the name of the input file from which the loop data were extracted.
#' }
#'
#' @importFrom utils read.table
#'
#' @examples
#' \dontrun{
#' loop_data <- loadloops(loop_file = c("loop1.bedpe", "loop2.bedpe"),
#' file_name = c("loop1", "loop2"),sep = "\t")
#' }
#'
#' @export
loadloops <- function(loop_file = NULL,
                      file_name = NULL,
                      sep = "\t"){
  plyr::ldply(1:length(loop_file),function(x){
    # read in
    bedpe <- read.table(loop_file[x],sep = sep,header = FALSE)

    # filename
    if(is.null(file_name)){
      fn <- loop_file[x]
    }else{
      fn <- file_name[x]
    }

    # check length
    if(ncol(bedpe) >= 6){
      bedpe <- bedpe[,1:6]
      colnames(bedpe) <- c("chrom1","start1","end1","chrom2","start2","end2")

      # return
      bedpe_df <- data.frame(seqnames = bedpe$chrom1,
                             start = (bedpe$start1 + bedpe$end1)/2,
                             end = (bedpe$start2 + bedpe$end2)/2,
                             score = (bedpe$end2 - bedpe$start1)/10^5,
                             fileName = fn)
    }else if(ncol(bedpe) == 3){
      bedpe_df <- bedpe
      bedpe_df$fileName <- fn
    }else{
      message("Please supply three columns bed or bedpe format data.")
      # break()
    }
    return(bedpe_df)
  }) -> loop_data
  return(loop_data)
}


#' Load Junction Data
#'
#' This function loads junction data from either BAM files or tab-delimited text
#' files.
#'
#' @param data_path A character vector of file paths to the BAM files or
#' tab-delimited text files.
#' @param file_name A character vector of file names for the BAM files or
#' tab-delimited text files. Optional if file names can be extracted from data_path.
#'
#' @return A data frame containing the junction data with columns for chromosome,
#' start position, end position, score, and file name.
#'
#' @examples
#' \dontrun{
#' # Load junctions from BAM files
#' bam_paths <- c("path/to/file1.bam", "path/to/file2.bam")
#' junctions <- loadJunction(data_path = bam_paths)
#'
#' # Load junctions from text files
#' txt_paths <- c("path/to/file1.txt", "path/to/file2.txt")
#' txt_names <- c("file1", "file2")
#' junctions <- loadJunction(data_path = txt_paths, file_name = txt_names)
#' }
#'
#' @export
loadJunction <- function(data_path = NULL,
                         file_name = NULL){
  # check data type
  if(endsWith(data_path[1],".bam")){
    # extract junctions
    plyr::ldply(seq_along(data_path),
                .progress = "text",
                function(x){
                  # file name
                  if(is.null(file_name)){
                    file_name_tmp <- unlist(strsplit(data_path[x],split = "/|.bam"))
                    file_name_ex <- file_name_tmp[length(file_name_tmp)]
                  }else{
                    file_name_ex <- file_name[x]
                  }

                  jc <- megadepth::bam_to_junctions(data_path[x])
                  jc_df <- megadepth::read_junction_table(jc)

                  # filter
                  data.table::setDT(jc_df)

                  jc_df_sum <- jc_df[unique == 1]
                  jc_df_sum <- jc_df_sum[, .(chr, start, end)]
                  jc_df_sum <- jc_df_sum[, .N, by = .(chr, start, end)]
                  jc_df_sum <- jc_df_sum[,fileName := file_name_ex]

                  return(jc_df_sum)
                }) -> junc_res
    colnames(junc_res) <- c("seqnames","start","end","score","fileName")
  }else{
    plyr::ldply(seq_along(data_path),
                .progress = "text",
                function(x){
                  # file name
                  if(is.null(file_name)){
                    message("please supply your file name with file_name params!")
                  }else{
                    file_name_ex <- file_name[x]
                  }

                  # read data
                  jc_df_sum <- utils::read.delim(data_path[x],header = FALSE,sep = "\t")
                  jc_df_sum <- jc_df_sum[,1:4]
                  jc_df_sum$fileName <- file_name_ex

                  return(jc_df_sum)
                }) -> junc_res
    colnames(junc_res) <- c("seqnames","start","end","score","fileName")
  }
  return(junc_res)
}
