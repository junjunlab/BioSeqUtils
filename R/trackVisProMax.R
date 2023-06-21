globalVariables(c("Freq","dist", "element_line", "exon_len", "facetted_pos_scales", "gene_group",
                  "gene_group2", "geom_label", "ggplotGrob", "guide_legend", "label_pos", "labeller",
                  "margin", "sampleName", "sample_group", "sample_group2", "scale_color_gradientn",
                  "scale_color_manual", "scale_fill_gradientn", "scale_x_continuous","smin_new",
                  "scale_y_continuous", "score", "smax", "smax_label", "smax_new", "smin",
                  "theme_void", "track_type", "xend", "yend", "ymax", "ymin"))

#' Visualize genomic data using BioSeqUtils package
#'
#' @author JunZhang
#'
#' This function generates a multi-track visualization of different genomic data
#' types, including gene annotation, signal tracks, loops, Hi-C interactions, bed
#' files, and junction files. For more details and full documentation, please refer
#' to [BioSeqUtils-manual](https://junjunlab.github.io/BioSeqUtils-manual/).
#'
#' @references [https://junjunlab.github.io/BioSeqUtils-manual/](https://junjunlab.github.io/BioSeqUtils-manual/)
#' @references [https://github.com/junjunlab/BioSeqUtils](https://github.com/junjunlab/BioSeqUtils)
#'
#' @param Input_gtf A GTF file containing the gene annotation information
#' with data frame format. Must be specified if you want to add gene annotation tracks.
#' @param Input_gene A character vectors. Must be specified if you want to add
#' @param upstream_extend Extend bases upstream for gene transcription start site,
#' default 0, accepting one or more vectors.
#' @param downstream_extend Extend bases downstream for gene transcription end site,
#' default 0, accepting one or more vectors.
#' gene annotation tracks.
#' @param Input_bw A bigWig file from "loadBigWig" function containing the signal
#' information for signal tracks.
#' @param Input_loop A loop file from "loadloops" function containing the loop
#' information for loop tracks.
#' @param Loop_curve_geom The curve type for loop tracks, "geom_arch2"(default)
#' or "geom_arch".
#' @param Input_hic A Hi-C interaction file from "prepareHic" function containing
#' the Hi-C interaction information for heatmap tracks.
#' @param Input_bed A bed file from "loadBed" function containing the bed
#' information for bed tracks.
#' @param Input_junction A junction file from "loadJunction" function containing
#' the junction information for junction tracks.
#' @param query_region A named list specifying the region to be visualized.
#' It must contain three elements: query_chr, query_start, and query_end.
#' @param signal_layer_bw_params A list of parameters for configuring signal
#' tracks created by bigWig files. Passed by "ggplot2::geom_rect" function.
#' @param signal_layer_loop_params A list of parameters for configuring loop
#' tracks created by loop files. Passed by "ggbio::geom_arch" or "jjPlot::geom_arch2"
#' function.
#' @param signal_layer_heatmap_params A list of parameters for configuring
#' heatmap tracks created by Hi-C interaction files. Passed by
#' "ggplot2::geom_polygon" function.
#' @param peaks_layer_params A list of parameters for configuring bed tracks
#' created by bed files. Passed by "ggplot2::geom_rect" function.
#' @param junc_layer_combined A boolean indicating whether to combine splice
#' junctions from multiple samples into a single track. Passed by
#' "ggplot2::geom_rect" or "ggplot2::geom_polygon" function.
#' @param add_band_line A boolean indicating whether to add band-line for junction
#' curve.
#' @param band_width A numeric value specifying the width of band line.
#' @param signal_layer_junction_params A list of parameters for configuring
#' junction tracks created by junction files. Passed by
#' "ggbio::geom_arch" or "ggplot2::geom_polygon" function.
#' @param signal_layer_junction_label_params A list of parameters for configuring
#' labels for junction tracks. Passed by "ggplot2::geom_label" function.
#' @param reverse_y_vars Character vectors indicating whether to reverse
#' the y-axis for certain tracks.
#' @param gene_group_info A named list containing group information for genes.
#' @param gene_group_info_order Character vectors specifying the order of gene groups.
#' @param gene_group_info2 A named list containing additional group information for genes.
#' @param gene_group_info2_order Character vectors specifying the order of
#' additional gene groups.
#' @param sample_group_info A named list containing group information for samples.
#' @param sample_group_info_order Character vectors specifying the order of sample groups.
#' @param sample_group_info2 A named list containing additional group information for samples.
#' @param sample_group_info2_order Character vectors specifying the order of
#' additional sample groups.
#' @param gene_order Character vectors specifying the order of genes in the visualization.
#' @param sample_order Character vectors specifying the order of samples in the visualization.
#' @param draw_chromosome A boolean indicating whether to draw the chromosome
#' ideogram on the bottom of the plot.
#' @param draw_chromosome_params A list of parameters for configuring the chromosome ideogram.
#' Details see "drawChromosome" function.
#' @param trans_topN An integer specifying the number of top transcripts to be shown per gene.
#' @param collapse_trans A boolean indicating whether to collapse exons from the same transcript into a single block.
#' @param exon_width A numeric value specifying the width of each exon block.
#' @param peak_width A numeric value specifying the width of each peak block.
#' @param add_gene_label_layer A boolean indicating whether to add gene labels to the plot.
#' @param gene_label_shift_y A numeric value specifying the vertical shift of gene labels.
#' @param gene_label_params A list of parameters for configuring gene labels.
#' Passed by "ggplot2::geom_text" function.
#' @param arrow_rel_len_params_list A list of parameters for configuring arrow
#' relative length. Details see "createSegment" function.
#' @param trans_exon_arrow_params A list of parameters for configuring arrows
#' of a transcript. Passed by "ggplot2::geom_segment" function.
#' @param trans_exon_col_params A list of parameters for coloring exons. Supplying
#' with "fill" or "color" parameters.
#' @param gene_dist_mark_params A list of parameters for control the segment line
#' in the "trans" track. Passed by "ggplot2::geom_segment" function.
#' @param gene_dist_mark_text_params A list of parameters for controling the
#' textl label in the "trans" track. Passed by "ggplot2::geom_text" function.
#' @param signal_range_pos A two elements vector specifying the position of signal
#' range label in each track.
#' @param signal_range A named list specifying the range of signal values to be
#' displayed in signal tracks.
#' @param panel_size_setting A list of settings for each panel size. Details see
#' "ggh4x::force_panelsizes".
#' @param fixed_column_range A logical value indicating whether to use a fixed
#' column range.
#' @param signal_range_label_params A list of parameters for signal range label.
#' Passed by "ggplot2::geom_text" function.
#' @param higlight_region Null or a named list of named values representing the
#' highlighted region. Like "higlight_region <- list(Actb = list(start =
#' c(142904000),end = c(142904500)))".
#' @param higlight_col Null or a named list of values representing the color of
#' the highlighted region. Like "higlight_col <- list(Actb = c("grey50"))".
#' @param higlight_col_alpha A numeric value indicating the alpha level of the
#' highlighted region color.
#' @param background_color_region Null or a named list of named values representing
#' the background color of the region. Like "background_color_region <- list(Actb
#' = c("Input" = "grey90","IP-1"="yellow","IP-2"="purple"),
#' Myc = c("Input" = "red","IP-1"="blue","IP-2"="green"))"
#' "background_color_region <- list("chr5:65438603-65834288" = c("Input" = "grey90",
#' "IP-1"="yellow","IP-2"="purple"),"chr15:9263295-9785835" = c("Input" = "red",
#' "IP-1"="blue","IP-2"="green"))"
#' @param background_region_alpha A numeric value indicating the alpha level of
#' the background color of the region.
#' @param by_layer_x A logical value indicating whether to plot the strip by layer
#' in the x-axis.
#' @param by_layer_y A logical value indicating whether to plot the strip by layer
#' in the y-axis.
#' @param column_strip_setting_list A list of settings for column strips. Passed
#' with "ggh4x::elem_list_rect".
#' @param row_strip_setting_list A list of settings for row strips. Passed with
#' "ggh4x::elem_list_rect".
#' @param column_strip_text_setting_list A list of settings for column strip text.
#' Passed with "ggh4x::elem_list_text".
#' @param row_strip_text_setting_list A list of settings for row strip text.
#' Passed with "ggh4x::elem_list_text".
#' @param add_gene_region_label A logical value indicating whether to add genomic
#' region labels under the gene symbols.
#' @param base_size A numeric value indicating the base font size.
#' @param panel.spacing A numeric vector of values indicating the spacing between panels.
#' @param sample_fill_col Null or a numeric vector of values representing the fill
#' color of the samples.
#' @param loops_col Null or character vector of values representing the color of
#' the loops.
#' @param heatmap_fill_col Null or a character vector of values representing the
#' fill color of the heatmap.
#' @param peak_fill_col Null or a character vector of values representing the fill
#' color of the peaks.
#' @param trans_fill_col Null or a character vector of values representing the fill
#' color of the trans.
#' @param remove_chrom_panel_border A logical value indicating whether to remove
#' the panel border for the chromosome track.
#' @param remove_all_panel_border A logical value indicating whether to remove
#' all panel borders.
#' @param xlimit_range Null or a numeric vector of values representing the x-axis
#' limit range.
#' @param Intron_line_type Line type for intron regions, "line"(default) or "chevron".
#' @param show_y_ticks Whether show Y axis ticks instead of range label, default FALSE.
#'
#'
#' @return GGPLOT
#' @export
trackVisProMax <- function(Input_gtf = NULL,
                           Input_gene = NULL,
                           upstream_extend = 0,
                           downstream_extend = 0,
                           Input_bw = NULL,
                           Input_loop = NULL,
                           Loop_curve_geom = "geom_arch2",
                           Input_hic = NULL,
                           Input_bed = NULL,
                           Input_junction = NULL,
                           query_region = list(query_chr = NULL,
                                               query_start = NULL,
                                               query_end = NULL),
                           signal_layer_bw_params = list(),
                           signal_layer_loop_params = list(),
                           signal_layer_heatmap_params = list(),
                           peaks_layer_params = list(),
                           junc_layer_combined = FALSE,
                           add_band_line = FALSE,
                           band_width = 0.5,
                           signal_layer_junction_params = list(),
                           signal_layer_junction_label_params = list(),
                           reverse_y_vars = NULL,
                           gene_group_info = NULL,
                           gene_group_info_order = NULL,
                           gene_group_info2 = NULL,
                           gene_group_info2_order = NULL,
                           sample_group_info = NULL,
                           sample_group_info_order = NULL,
                           sample_group_info2 = NULL,
                           sample_group_info2_order = NULL,
                           gene_order = NULL,
                           sample_order = NULL,
                           draw_chromosome = FALSE,
                           draw_chromosome_params = list(ideogram_obj = NULL),
                           trans_topN = 2,
                           collapse_trans = FALSE,
                           exon_width = 0.5,
                           peak_width = 0.5,
                           add_gene_label_layer = FALSE,
                           gene_label_shift_y = 0.3,
                           gene_label_params = list(),
                           arrow_rel_len_params_list = list(),
                           trans_exon_arrow_params = list(),
                           trans_exon_col_params = list(),
                           gene_dist_mark_params = list(),
                           gene_dist_mark_text_params = list(),
                           signal_range_pos = c(0.85,0.85),
                           signal_range = NULL,
                           panel_size_setting = list(),
                           fixed_column_range = TRUE,
                           signal_range_label_params = list(),
                           higlight_region = NULL,
                           higlight_col = NULL,
                           higlight_col_alpha = 0.2,
                           background_color_region = NULL,
                           background_region_alpha = 0.2,
                           by_layer_x = FALSE,
                           by_layer_y = FALSE,
                           column_strip_setting_list = list(),
                           row_strip_setting_list = list(),
                           column_strip_text_setting_list = list(),
                           row_strip_text_setting_list = list(),
                           add_gene_region_label = FALSE,
                           base_size = 14,
                           panel.spacing = c(0.2,0),
                           sample_fill_col = NULL,
                           loops_col = NULL,
                           heatmap_fill_col = NULL,
                           peak_fill_col = NULL,
                           trans_fill_col = NULL,
                           remove_chrom_panel_border = FALSE,
                           remove_all_panel_border = FALSE,
                           xlimit_range = NULL,
                           Intron_line_type = "line",
                           show_y_ticks = FALSE
){
  options(warn=-1)
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  # ==============================================================================
  # 1_dplyr::filter target gene region for signals
  # ==============================================================================
  gtf <- Input_gtf
  bw <- Input_bw
  bed <- Input_bed
  loop <- Input_loop
  heatmap <- Input_hic
  junction <- Input_junction

  # Input_bed <- bed_df
  # bed <- Input_bed

  if(!is.null(Input_bw)){
    bw$track_type <- "bigwig"
    bw$id <- NA
  }

  if(!is.null(Input_loop)){
    loop$track_type <- "loop"
    loop$id <- NA
  }

  if(!is.null(Input_hic)){
    heatmap$track_type <- "heatmap"
  }

  if(!is.null(Input_junction)){
    junction$track_type <- "junction"
    junction$id <- NA
  }

  if(!is.null(Input_bw) | !is.null(Input_loop) | !is.null(Input_hic) | !is.null(Input_junction)){
    # input_signal_file <- rbind(bw,loop,heatmap,junction)
    input_signal_file <- plyr::rbind.fill(bw,loop,heatmap,junction)
  }else{
    input_signal_file <- data.frame(seqnames = NA,start = NA,end = NA,score = NA,
                                    fileName = "trans",track_type = NA,id = NA)
  }

  # =====================================================
  # dplyr::filter signal data
  # x = 1
  if(!is.null(Input_bw)){
    if(is.null(Input_gene)){
      plyr::ldply(1:length(query_region[[1]]),function(x){
        # dplyr::filter signals
        tmp <- input_signal_file %>%
          dplyr::filter(track_type != "heatmap") %>%
          dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
          dplyr::filter(start >= query_region$query_start[x] & end <= query_region$query_end[x])

        tmp_cont <- input_signal_file %>%
          dplyr::filter(track_type == "heatmap") %>%
          dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
          dplyr::filter(start >= query_region$query_start[x] & start <= query_region$query_end[x])

        # remove no 4 points position
        id_n <- data.frame(table(tmp_cont$id)) %>% dplyr::filter(Freq == 4)
        tmp_cont <- tmp_cont[which(tmp_cont$id %in% id_n$Var1),]

        mer_tmp <- plyr::rbind.fill(tmp,tmp_cont) %>%
          mutate(gene = paste(ifelse(startsWith(query_region$query_chr[x],"chr"),
                                     query_region$query_chr[x],
                                     paste("chr",query_region$query_chr[x],sep = "")),
                              ":",as.integer(query_region$query_start[x]),"-",
                              as.integer(query_region$query_end[x]),sep = ""))
      }) -> region.df

    }else(
      # x = 1
      # upstream_extend = 0
      # downstream_extend = 0
      plyr::ldply(seq_along(Input_gene),function(x){
        tmp <- gtf %>%
          dplyr::filter(gene_name == Input_gene[x])
        chr <- as.character(unique(tmp$seqnames))
        xmin = min(tmp$start) - ifelse(length(upstream_extend) == 1,upstream_extend,upstream_extend[x])
        xmax = max(tmp$end) + ifelse(length(downstream_extend) == 1,downstream_extend,downstream_extend[x])

        # dplyr::filter signals
        sig <- input_signal_file %>%
          dplyr::filter(seqnames %in% chr) %>%
          dplyr::filter(start >= xmin & end <= xmax) %>%
          mutate(gene = Input_gene[x])
        return(sig)
      }) -> region.df
    )
  }else{
    # region.df <- NULL
    if(is.null(Input_gene)){
      plyr::ldply(1:length(query_region[[1]]),function(x){
        # dplyr::filter signals
        if(is.null(Input_loop) & is.null(Input_hic)){
          tmp_tmp <- input_signal_file
        }else if(!is.null(Input_loop) & is.null(Input_hic)){
          tmp <- input_signal_file %>%
            dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
            dplyr::filter(start >= query_region$query_start[x] & end <= query_region$query_end[x])
          tmp_tmp <- tmp
        }else if(!is.null(Input_hic)){
          tmp <- input_signal_file %>%
            dplyr::filter(track_type != "heatmap") %>%
            dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
            dplyr::filter(start >= query_region$query_start[x] & end <= query_region$query_end[x])

          tmp_cont <- input_signal_file %>%
            dplyr::filter(track_type == "heatmap") %>%
            dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
            dplyr::filter(start >= query_region$query_start[x] & start <= query_region$query_end[x])

          # remove no 4 points position
          id_n <- data.frame(table(tmp_cont$fileName,tmp_cont$id)) %>% dplyr::filter(Freq == 4)
          tmp_cont <- tmp_cont[which(tmp_cont$id %in% unique(id_n$Var2)),]

          tmp_tmp <- plyr::rbind.fill(tmp,tmp_cont)
        }

        mer_tmp <- tmp_tmp %>%
          mutate(gene = paste(ifelse(startsWith(query_region$query_chr[x],"chr"),
                                     query_region$query_chr[x],
                                     paste("chr",query_region$query_chr[x],sep = "")),
                              ":",as.integer(query_region$query_start[x]),"-",
                              as.integer(query_region$query_end[x]),sep = ""))
      }) -> region.df
    }else{
      plyr::ldply(seq_along(Input_gene),function(x){
        # dplyr::filter signals
        sig <- input_signal_file %>%
          mutate(gene = Input_gene[x])
        return(sig)
      }) -> region.df
    }
  }

  # =====================================================
  # dplyr::filter peaks files
  if(!is.null(Input_gene) & !is.null(Input_bed)){
    plyr::ldply(seq_along(Input_gene),function(x){
      tmp <- gtf %>%
        dplyr::filter(gene_name == Input_gene[x])
      chr <- as.character(unique(tmp$seqnames))
      xmin = min(tmp$start) - ifelse(length(upstream_extend) == 1,upstream_extend,upstream_extend[x])
      xmax = max(tmp$end) + ifelse(length(downstream_extend) == 1,downstream_extend,downstream_extend[x])

      # dplyr::filter signals
      # peak_width = 0.5
      sig <- bed %>%
        dplyr::filter(seqnames %in% chr) %>%
        dplyr::filter(start >= xmin & end <= xmax) %>%
        mutate(gene = Input_gene[x],
               ymin = y - peak_width*0.5,
               ymax = y + peak_width*0.5)
      return(sig)
    }) -> bed.df
  }else if(is.null(Input_gene) & !is.null(Input_bed)){
    plyr::ldply(1:length(query_region[[1]]),function(x){
      # dplyr::filter signals
      tmp <- bed %>%
        dplyr::filter(seqnames %in% query_region$query_chr[x]) %>%
        dplyr::filter(start >= query_region$query_start[x] & end <= query_region$query_end[x]) %>%
        mutate(gene = paste(ifelse(startsWith(query_region$query_chr[x],"chr"),
                                   query_region$query_chr[x],
                                   paste("chr",query_region$query_chr[x],sep = "")),
                            ":",query_region$query_start[x],"-",query_region$query_end[x],sep = ""),
               ymin = y - peak_width*y*0.5,
               ymax = y + peak_width*y*0.5)
    }) -> bed.df
  }

  # ==============================================================================
  # 2_add trans and chromo facet
  # ==============================================================================
  # draw_chromosome = F
  if(!is.null(Input_bed)){
    add_facet_name = c("peaks","trans")
  }else{
    add_facet_name = c("trans")
  }

  if(draw_chromosome == TRUE){
    add_facet_name = append(add_facet_name,c("chrom"))
  }

  # add info data
  if(is.null(Input_gene)){
    gene = unique(region.df$gene)
  }else{
    gene = Input_gene
  }
  tran_facet <- data.frame(seqnames = NA,start = NA,end = NA,score = NA,
                           fileName = rep(add_facet_name,each = length(unique(region.df$gene))),
                           track_type = NA,
                           id = NA,
                           gene = gene)

  # merge
  mer.df <- plyr::rbind.fill(region.df,tran_facet) %>% unique()

  # ==============================================================================
  # 3_add group info for genes or samples
  # ==============================================================================
  gene_group_info = gene_group_info
  gene_group_info2 = gene_group_info2
  sample_group_info = sample_group_info
  sample_group_info2 = sample_group_info2

  # sample_group_info = list("control" = "Input",
  #                          "treat" = c("IP-1","IP-2"))
  # sample_group_info2 = list("sample1" = c("Input","IP-1"),
  #                          "sample2" = c("IP-2"))
  # gene_group_info = NULL
  # sample_group_info = NULL

  # add gene group
  if(!is.null(gene_group_info)){
    # g = 1
    # add gene group
    plyr::ldply(1:length(gene_group_info),function(g){
      tmp <- mer.df[which(mer.df$gene %in% gene_group_info[[g]]),]
      tmp <- tmp %>%
        mutate(gene_group = names(gene_group_info[g]))
    }) -> tmp1

    # orders
    if(is.null(gene_group_info_order)){
      gene_levels <- as.character(names(gene_group_info))
    }else{
      gene_levels <- gene_group_info_order
    }

    tmp1$gene_group <- factor(tmp1$gene_group,levels = gene_levels)
  }else{
    tmp1 <- mer.df %>%
      mutate(gene_group = NA)
  }

  # tmp1 %>% dplyr::select(gene,gene_group) %>% distinct()

  if(!is.null(gene_group_info2)){
    # g = 1
    # add gene group
    plyr::ldply(1:length(gene_group_info2),function(g){
      tmp <- tmp1[which(tmp1$gene %in% gene_group_info2[[g]]),]
      tmp <- tmp %>%
        mutate(gene_group2 = names(gene_group_info2[g]))
    }) -> tmp1

    # tmp1 %>% dplyr::select(gene,gene_group,gene_group2) %>% distinct()

    # orders
    if(is.null(gene_group_info2_order)){
      gene_levels <- as.character(names(gene_group_info2))
    }else{
      gene_levels <- gene_group_info2_order
    }

    tmp1$gene_group2 <- factor(tmp1$gene_group2,levels = gene_levels)
  }else{
    tmp1 <- tmp1 %>%
      mutate(gene_group2 = NA)
  }

  # ================================================================
  # add sample group
  if(!is.null(sample_group_info)){
    sample_group_info = append(sample_group_info,list(other = add_facet_name))
    # s = 1
    # add sample group
    plyr::ldply(1:length(sample_group_info),function(s){
      tmp <- tmp1[which(tmp1$fileName %in% sample_group_info[[s]]),]
      tmp <- tmp %>%
        mutate(sample_group = names(sample_group_info[s]))
    }) -> tmp2

    # orders
    if(is.null(sample_group_info_order)){
      sample_levels <- as.character(names(sample_group_info))
    }else{
      sample_levels <- sample_group_info_order
    }

    tmp2$sample_group <- factor(tmp2$sample_group,levels = sample_levels)
  }else{
    tmp2 <- tmp1 %>%
      mutate(sample_group = NA)
  }

  if(!is.null(sample_group_info2)){
    sample_group_info2 = append(sample_group_info2,list(other = add_facet_name))
    # s = 1
    # add sample group
    plyr::ldply(1:length(sample_group_info2),function(s){
      tmp <- tmp2[which(tmp2$fileName %in% sample_group_info2[[s]]),]
      tmp <- tmp %>%
        mutate(sample_group2 = names(sample_group_info2[s]))
    }) -> tmp2

    # orders
    if(is.null(sample_group_info2_order)){
      sample_levels <- as.character(names(sample_group_info2))
    }else{
      sample_levels <- sample_group_info2_order
    }

    tmp2$sample_group2 <- factor(tmp2$sample_group2,levels = sample_levels)
  }else{
    tmp2 <- tmp2 %>%
      mutate(sample_group2 = NA)
  }

  # gene orders
  # gene_order = NULL
  if(is.null(gene_order)){
    tmp2$gene <- factor(tmp2$gene,levels = unique(region.df$gene))
  }else{
    tmp2$gene <- factor(tmp2$gene,levels = gene_order)
  }

  # sample orders
  # sample_order = NULL
  if(is.null(sample_order)){
    tmp2$fileName <- factor(tmp2$fileName,
                            levels = unique(c(unique(region.df$fileName),add_facet_name)))
  }else{
    tmp2$fileName <- factor(tmp2$fileName,levels = c(sample_order,add_facet_name))
  }

  # str(tmp2)

  # x = 1
  if(!is.null(Input_bed)){
    plyr::ldply(1:nrow(bed.df),function(x){
      tmp <- bed.df[x,]
      tmp <- tmp %>%
        mutate(fileName = factor("peaks",levels = levels(tmp2$fileName)),
               gene = factor(gene,levels = levels(tmp2$gene)),
               gene_group = factor(unique(tmp2[which(tmp2$gene == gene),"gene_group"]),
                                   levels = levels(tmp2$gene_group)),
               gene_group2 = factor(unique(tmp2[which(tmp2$gene == gene),"gene_group2"]),
                                    levels = levels(tmp2$gene_group2)),
               sample_group = factor(unique(tmp2[which(tmp2$gene == gene),"sample_group"]),
                                     levels = levels(tmp2$sample_group)),
               sample_group2 = factor(unique(tmp2[which(tmp2$gene == gene),"sample_group2"]),
                                      levels = levels(tmp2$sample_group2)))

    }) -> bed.df.new
  }

  # str(bed.df.new)
  # tmp2 %>% dplyr::select(gene,gene_group,gene_group2,
  #                 sample_group,sample_group2) %>% distinct()

  # ==============================================================================
  # 4_gene structures
  # ==============================================================================
  # tmp2 %>% dplyr::select(gene,gene_group,gene_group2) %>% distinct() %>%
  #   arrange(gene)
  # x = 1
  # collapse_trans = FALSE
  # x = 1
  if(is.null(Input_gene)){
    tmp_gtf <- plyr::ldply(1:length(query_region[[1]]),function(x){
      tmp <- gtf %>%
        dplyr::filter(seqnames %in% query_region$query_chr[x] &
                        start >= query_region$query_start[x] &
                        end <= query_region$query_end[x] &
                        type != "gene") %>%
        mutate(gene = paste(ifelse(startsWith(query_region$query_chr[x],"chr"),
                                   query_region$query_chr[x],
                                   paste("chr",query_region$query_chr[x],sep = "")),
                            ":",as.integer(query_region$query_start[x]),"-",
                            as.integer(query_region$query_end[x]),sep = ""))
      return(tmp)
    })
  }else{
    tmp_gtf <- gtf %>%
      dplyr::filter(gene_name %in% Input_gene & type != "gene") %>%
      mutate(gene = gene_name)
  }

  # get gene id
  gid <- unique(tmp_gtf$gene_id)

  # x = 1
  # expand gene
  purrr::map_df(1:length(gid),function(x){
    tmp <- tmp_gtf %>% dplyr::filter(gene_id == gid[x])

    # whether has transcript in gene info
    if('transcript' %in% unique(tmp$type)){
      new_tmp <- tmp
    }else{
      # add transcript feature
      tid_n <- unique(tmp$transcript_id)
      new <- tmp %>% group_by(seqnames,gene,gene_id,gene_name,transcript_id,strand) %>%
        summarise(start = min(start),end = max(end)) %>%
        mutate(type = "transcript")
      new_tmp <- tmp %>% add_row(seqnames = new$seqnames,
                                 start = new$start,end = new$end,
                                 type = new$type,
                                 transcript_id = new$transcript_id,
                                 gene = new$gene,
                                 gene_id = new$gene_id,
                                 gene_name = new$gene_name,
                                 strand = new$strand)
    }

    # dplyr::filter top transcripts
    if("exon" %in% new_tmp$type){
      get_type = "exon"
    }else{
      get_type = "CDS"
    }

    trans_len <- new_tmp %>%
      dplyr::filter(type == get_type) %>%
      dplyr::group_by(transcript_id,type) %>%
      dplyr::summarise(exon_len = sum(width)) %>%
      dplyr::arrange(desc(exon_len)) %>% ungroup()

    # choose numbers transcript to show
    if(trans_topN == "all"){
      trans_len <- trans_len
    }else{
      trans_len <- trans_len %>% slice_head(.,n = trans_topN)
    }

    filtered_trans <- new_tmp %>%
      # dplyr::filter(transcript_id %in% trans_len$transcript_id & type != "transcript")
      dplyr::filter(transcript_id %in% trans_len$transcript_id)

    # x = 1
    tid = trans_len$transcript_id
    plyr::ldply(seq_along(tid),function(x){
      # whether collapse gene
      if(collapse_trans == TRUE){
        y_p = 1
        trans_topN <<- 1
      }else{
        y_p = x
        # if(!is.null(Input_gene)){
        #   trans_topN <- nrow(trans_len)
        # }else{
        #   trans_topN <- trans_topN
        # }
      }

      tmp <- filtered_trans %>% dplyr::filter(transcript_id %in% rev(tid)[x]) %>%
        mutate(ymin = if_else(type %in% c("5UTR","five_prime_utr","3UTR","three_prime_utr"),
                              y_p - exon_width*0.25,y_p - exon_width*0.5),
               ymax = if_else(type %in% c("5UTR","five_prime_utr","3UTR","three_prime_utr"),
                              y_p + exon_width*0.25,y_p + exon_width*0.5),
               y = y_p)

      # remove exon if it is a coding gene
      if("CDS" %in% unique(tmp$type)){
        tmp <- tmp %>% dplyr::filter(type != "exon")
      }else{
        tmp <- tmp
      }
      return(tmp)
    }) -> trans_pos

    # add gene_group_info sample_group_info
    if(is.null(Input_gene)){
      # idx <- match(unique(trans_pos$seqnames),query_region$query_chr)
      # trans_pos <- trans_pos %>% mutate(gene = levels(tmp2$gene)[idx],fileName = "trans")
      trans_pos <- trans_pos %>% mutate(fileName = "trans")
    }else{
      trans_pos <- trans_pos %>% mutate(gene = gene_name,fileName = "trans")
    }

    if(!is.null(sample_group_info)){
      trans_pos$sample_group <- "other"

      # orders
      sample_levels <- as.character(names(sample_group_info))
      trans_pos$sample_group <- factor(trans_pos$sample_group,levels = sample_levels)
    }else{
      trans_pos$sample_group <- NA
    }

    if(!is.null(sample_group_info2)){
      trans_pos$sample_group2 <- "other"

      # orders
      sample_levels <- as.character(names(sample_group_info2))
      trans_pos$sample_group2 <- factor(trans_pos$sample_group2,levels = sample_levels)
    }else{
      trans_pos$sample_group2 <- NA
    }

    if(!is.null(gene_group_info)){
      gene.group.tmp <- tmp2[,c("gene","gene_group")] %>% unique()
      trans_pos$gene_group <- gene.group.tmp[which(gene.group.tmp$gene == unique(trans_pos$gene)),
                                             "gene_group"]

      # orders
      gene_levels <- as.character(names(gene_group_info))
      trans_pos$gene_group <- factor(trans_pos$gene_group,levels = gene_levels)
    }else{
      trans_pos$gene_group <- NA
    }

    if(!is.null(gene_group_info2)){
      gene.group.tmp <- tmp2[,c("gene","gene_group2")] %>% unique()

      trans_pos <- trans_pos %>% left_join(.,y = gene.group.tmp,by = "gene")
      # trans_pos$gene_group2 <- gene.group.tmp[which(gene.group.tmp$gene == unique(trans_pos$gene)),
      #                                         "gene_group2"]

      # orders
      gene_levels <- as.character(names(gene_group_info2))
      trans_pos$gene_group2 <- factor(trans_pos$gene_group2,levels = gene_levels)
    }else{
      trans_pos$gene_group2 <- NA
    }

    return(trans_pos)
    # return(exon_ypos)
  }) -> transcript.df

  # order
  transcript.df$fileName <- factor(transcript.df$fileName,levels = levels(tmp2$fileName))
  transcript.df$gene <- factor(transcript.df$gene,levels = levels(tmp2$gene))

  # =============================================
  # gene label
  # =============================================
  # add_gene_label_layer = FALSE
  if(add_gene_label_layer == TRUE){
    # gene_label_df
    gene_label_df <- transcript.df %>%
      dplyr::filter(type == "transcript")

    # gene label layer
    # gene_label_params = list(size = 3,label_aes = "gene_name",color = "black")
    # gene_label_layer <- geom_text_repel(data = gene_label_df,
    #                                     aes(x = (start + end)/2,y = y - exon_width/2 + 0.3,
    #                                         label = get(gene_label_params$label_aes)),
    #                                     min.segment.length = 0,
    #                                     max.overlaps = Inf,
    #                                     hjust = 0.5,
    #                                     size = gene_label_params$size,
    #                                     color = gene_label_params$color)

    if(collapse_trans == TRUE){
      gene_label_df <- gene_label_df %>%
        group_by(seqnames,gene_name,gene,y,fileName,sample_group,sample_group2,
                 gene_group,gene_group2) %>%
        summarise(start = min(start),end = max(end))
    }

    mapping = aes(x = (start + end)/2,y = y - exon_width/2 + gene_label_shift_y,
                  label = gene_name)

    # dplyr::select function
    if(!is.null(Input_gene)){
      gene_label_layer <- do.call(geom_text,
                                  modifyList(list(data = gene_label_df,
                                                  mapping = mapping,
                                                  hjust = 0.5,
                                                  size = 3,
                                                  color = "black"),gene_label_params))
    }else{
      gene_label_layer <- do.call(ggrepel::geom_text_repel,
                                  modifyList(list(data = gene_label_df,
                                                  mapping = mapping,
                                                  min.segment.length = 0,
                                                  max.overlaps = Inf,
                                                  hjust = 0.5,
                                                  size = 3,
                                                  color = "black"),gene_label_params))
    }
  }else{
    gene_label_layer <- NULL
  }

  # =============================================
  # transcript arrow
  # =============================================
  # x = 1
  plyr::ldply(seq_along(gid),function(x){
    tmp <- transcript.df %>%
      dplyr::filter(gene_id == gid[x] & type == "transcript")

    # create segment
    plyr::ldply(1:nrow(tmp),function(x){
      tmp1 <- tmp[x,]
      ypos = unique(transcript.df[which(transcript.df$transcript_id == tmp1$transcript_id),
                                  c("y")])
      # arrow_rel_len_params_list
      # arrow_rel_len_params_list <- list(n_division = NULL,rel_len = 0.1,abs_len = NULL)

      # choose intron line type
      # Intron_line_type = "chevron"
      if(Intron_line_type == "chevron"){
        tmp_exon <- gtf %>%
          filter(transcript_id ==  tmp1$transcript_id & type == "exon") %>%
          arrange(start,end)

        xstart = tmp_exon$end[1:nrow(tmp_exon) - 1]
        xend = tmp_exon$start[2:nrow(tmp_exon)]
        seg_data <- data.frame(x = c(xstart,(xstart + xend)/2),
                               xend = c((xstart + xend)/2,xend),
                               y = rep(c(ypos,ypos + exon_width*0.25),each = nrow(tmp_exon) - 1),
                               yend = rev(rep(c(ypos,ypos + exon_width*0.25),each = nrow(tmp_exon) - 1)),
                               transcript_id = tmp1$transcript_id
        )
      }else if(Intron_line_type == "line"){
        seg_data <-
          do.call(createSegment,modifyList(list(xPos = c(tmp1$start,tmp1$end),
                                                yPos = rep(ypos,2),
                                                rel_len = 0.08),
                                           arrow_rel_len_params_list)) %>%
          mutate(transcript_id = tmp1$transcript_id)
      }

      # get group info
      seg_data$gene <- tmp1$gene
      seg_data$ends <- ifelse(unique(tmp$strand) == "+","last","first")
      seg_data$fileName <- "trans"

      seg_data <- seg_data %>% mutate(gene_group = unique(tmp1$gene_group),
                                      gene_group2 = unique(tmp1$gene_group2),
                                      sample_group = unique(tmp1$sample_group),
                                      sample_group2 = unique(tmp1$sample_group2))

      return(seg_data)
    }) -> seg_arrow

    return(seg_arrow)
  }) -> final_arrow_data

  # order
  final_arrow_data$fileName <- factor(final_arrow_data$fileName,levels = levels(tmp2$fileName))
  final_arrow_data$gene <- factor(final_arrow_data$gene,levels = levels(tmp2$gene))

  # arrow layer stands strand information
  # x = "Myc"
  # trans_exon_arrow_params = list(length = 1,fill = "grey60",color = "grey60",linewidth = 0.75)
  trans_arrow_layer <- lapply(unique(final_arrow_data$transcript_id), function(x){
    tmp <- final_arrow_data[which(final_arrow_data$transcript_id == x),]

    # arrow layer
    if(Intron_line_type == "chevron"){
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = x,xend = xend,
                                  y = y,yend = yend),
                              linewidth = 0.75,
                              color = "grey60"),trans_exon_arrow_params))
    }else{
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = x,xend = xend,
                                  y = y,yend = yend),
                              linewidth = 0.75,
                              arrow = arrow(type = "closed",
                                            length = unit(1,"mm"),
                                            ends = unique(tmp$ends)),
                              arrow.fill = "grey60",
                              color = "grey60"),trans_exon_arrow_params))
    }
  })

  # trans_struct_layer
  # trans_exon_col_params = list(fill = "orange",color = "grey60")
  # trans_struct_layer <- geom_rect(data = transcript.df %>% dplyr::filter(type != "transcript"),
  #                                 aes(xmin = start,xmax = end,
  #                                     ymin = ymin, ymax = ymax),
  #                                 fill = trans_exon_col_params$fill,
  #                                 color = trans_exon_col_params$color)

  if(!is.null(trans_exon_col_params$mapping)){
    trans_mapping <- list(data = transcript.df %>% dplyr::filter(type != "transcript"),
                          mapping= aes(xmin = start,xmax = end,
                                       ymin = ymin, ymax = ymax,
                                       fill = "orange"),
                          color = "grey60")
  }else{
    trans_mapping <- list(data = transcript.df %>% dplyr::filter(type != "transcript"),
                          mapping= aes(xmin = start,xmax = end,
                                       ymin = ymin, ymax = ymax),
                          fill = "orange",
                          color = "grey60")
  }



  trans_struct_layer <- do.call(geom_rect,
                                modifyList(trans_mapping,trans_exon_col_params))

  # ==============================================================================
  # 5_segment and arrow data for chromosome label and region length
  # ==============================================================================
  if(!is.null(Input_gene)){
    segment.df <- transcript.df %>%
      group_by(fileName,gene,seqnames,strand,sample_group,sample_group2,gene_group,gene_group2)
  }else{
    segment.df <- transcript.df %>%
      group_by(fileName,gene,seqnames,sample_group,sample_group2,gene_group,gene_group2)
  }

  # add segment info
  segment.df <- segment.df %>%
    summarise(start = min(start),end = max(end))

  # whether re-assign x limits
  if(!is.null(xlimit_range)){
    if(length(xlimit_range) == 2 & !is.list(xlimit_range)){
      segment.df <- segment.df %>%
        mutate(start = xlimit_range[1],end = xlimit_range[2])
    }else{
      segment.df <- plyr::ldply(1:nrow(segment.df),function(x){
        tmp <- segment.df[x,]
        tmp <- tmp %>%
          mutate(start = xlimit_range[[x]][1],end = xlimit_range[[x]][2])
        return(tmp)
      })
    }
  }

  # whether re-assign x limits for upstream_extend and downstream_extend
  if(length(upstream_extend) == 1 & length(downstream_extend) == 1){
    segment.df <- segment.df %>%
      mutate(start = start - upstream_extend,end = end + downstream_extend)
  }else if(length(upstream_extend) == 1 & length(downstream_extend) >= 1){
    segment.df <- plyr::ldply(1:nrow(segment.df),function(x){
      tmp <- segment.df[x,]
      tmp <- tmp %>%
        mutate(start = start - upstream_extend,end = end + downstream_extend[x])
      return(tmp)
    })
  }else if(length(upstream_extend) >= 1 & length(downstream_extend) == 1){
    segment.df <- plyr::ldply(1:nrow(segment.df),function(x){
      tmp <- segment.df[x,]
      tmp <- tmp %>%
        mutate(start = start - upstream_extend[x],end = end + downstream_extend)
      return(tmp)
    })
  }else{
    segment.df <- plyr::ldply(1:nrow(segment.df),function(x){
      tmp <- segment.df[x,]
      tmp <- tmp %>%
        mutate(start = start - upstream_extend[x],end = end + downstream_extend[x])
      return(tmp)
    })
  }

  # two segment lines position
  segment.df <- segment.df %>%
    mutate(ar1_end = start + (end - start)/3,
           ar2_start = end - (end - start)/3)

  # arrow data on gene structures
  # x = 1
  plyr::ldply(1:nrow(segment.df),function(x){
    tmp <- segment.df[x,] %>%
      # add chr prefix
      mutate(seqnames = if_else(startsWith(as.character(seqnames),"chr"),
                                seqnames,paste("chr",seqnames,sep = "")))
    # calculate arrow positions
    # t_num = table(data.frame(unique(transcript.df[,c("gene","transcript_id")]))$gene)
    t_num = table(data.frame(unique(final_arrow_data[,c("gene","y")]))$gene)

    res <- data.frame(start = c(tmp$start,tmp$end - (tmp$end - tmp$start)/3),
                      end = c(tmp$start + (tmp$end - tmp$start)/3,tmp$end)) %>%
      mutate(fileName = tmp$fileName,gene = tmp$gene,strand = tmp$strand,
             .before = start) %>%
      mutate(label = paste(tmp$seqnames,": ",round((tmp$end - tmp$start)/10^3,digits = 2),
                           " kb",sep = ""),
             label_pos = (tmp$end + tmp$start)/2) %>%
      mutate(gene_group = tmp$gene_group,
             gene_group2 = tmp$gene_group2,
             sample_group = tmp$sample_group,
             sample_group2 = tmp$sample_group2,
             y = t_num[tmp$gene] + 1)

    # add arrow direction
    if(!is.null(Input_gene)){
      res <- res %>%
        mutate(ends = if_else(strand == "+","last","first"))
    }else{
      res <- res %>%
        mutate(ends = c("first","last"))
    }

    return(res)
  }) -> arrow.df

  # arrow layer stands strand information
  # gene_dist_mark_params = list(linewidth = 0.75,color = "black",length = 2)
  seg_arrow_layer <- lapply(unique(arrow.df$gene), function(x){
    tmp <- arrow.df %>% dplyr::filter(gene == x)
    # geom_segment(data = tmp,
    #              aes(x = start,xend = end,
    #                  y = y,yend = y),linewidth = gene_dist_mark_params$linewidth,
    #              color = gene_dist_mark_params$color,
    #              arrow = arrow(ends = tmp$ends,
    #                            length = unit(gene_dist_mark_params$length,"mm")))
    seg_ar <-
      do.call(geom_segment,
              modifyList(list(data = tmp,
                              aes(x = start,xend = end,
                                  y = if(collapse_trans == TRUE){y = 2}else{y},
                                  yend = if(collapse_trans == TRUE){y = 2}else{y}
                              ),
                              linewidth = 0.3,
                              color = "black",
                              arrow = arrow(ends = tmp$ends,
                                            length = unit(2,"mm"))),
                         gene_dist_mark_params))
    return(seg_ar)
  })

  # text label layer
  # gene_dist_mark_text_params = list(color = "black",size = 4)
  # text_layer <- geom_text(data = arrow.df,
  #                         aes(x = label_pos,y = y,label = label),
  #                         color = gene_dist_mark_text_params$color,
  #                         size = gene_dist_mark_text_params$size)

  text_layer <- do.call(geom_text,
                        modifyList(list(data = arrow.df,
                                        aes(x = label_pos,
                                            y = if(collapse_trans == TRUE){y = 2}else{y},
                                            label = label),
                                        color = "black",
                                        size = 3),
                                   gene_dist_mark_text_params))

  # # return layer
  # if(!is.null(Input_gene)){
  #
  # }else{
  #   seg_arrow_layer <- NULL
  #   text_layer <- NULL
  # }

  # ==============================================================================
  # 6_text and signal ranges for each panel
  # ==============================================================================
  if(!is.null(Input_bw) | !is.null(Input_loop) |!is.null(Input_hic) | !is.null(Input_junction)){
    # signal_range_pos = c(0.9,0.9)
    rg_xpos <- tmp2 %>%
      dplyr::filter(!(fileName %in% add_facet_name)) %>%
      group_by(gene) %>%
      summarise(xmin = min(start),xmax = max(end)) %>%
      ungroup() %>% mutate(rg_xpos = (xmax - xmin)*signal_range_pos[1] + xmin) %>%
      dplyr::select(gene,rg_xpos)

    # primitive range info for heatmap
    hetamp_y <- tmp2 %>% dplyr::filter(track_type == "heatmap") %>%
      group_by(fileName) %>% summarise(smax = max(end))

    # merge with rg_xpos
    rg_info <- tmp2 %>%
      dplyr::filter(!(fileName %in% add_facet_name)) %>%
      group_by(fileName,gene,track_type,gene_group,gene_group2,sample_group,sample_group2) %>%
      summarise(smin = min(score),smax = max(score)) %>%
      ungroup() %>%
      # mutate(rg_ypos = (smax - smin)*signal_range_pos[2] + smin) %>%
      left_join(.,rg_xpos,by = "gene")

    # whether dplyr::filter junction
    if(junc_layer_combined == TRUE){
      rg_info <- rg_info[which(rg_info$track_type != "junction"),]
    }

    # add smax for heatmap
    if(!is.null(Input_hic)){
      rg_info <- plyr::ldply(1:nrow(rg_info),function(x){
        tmp <- rg_info[x,]
        if(tmp$track_type == "heatmap"){
          tmp$smax <- hetamp_y[which(hetamp_y$fileName == tmp$fileName),]$smax
        }
        return(tmp)
      })
    }

    # whether supply new ranges
    # signal_range = list(Actb = c("Input" = 800,"IP-1" = 800,"IP-2" = 800),
    #                     Myc = c("Input" = 3,"IP-1" = 3,"IP-2" = 3))
    # signal_range_pos = c(0.9,0.9)
    # signal_range = NULL
    # fixed_column_range = TRUE

    # add new range column
    # x = 1
    if(!is.null(signal_range)){
      plyr::ldply(1:length(signal_range),function(x){
        # range new
        if(is.list(signal_range[[x]])){
          rg_tmp <-
            signal_range[[x]] %>%
            data.frame(check.names = FALSE) %>%
            t() %>% data.frame(check.names = FALSE) %>%
            tibble::rownames_to_column(var = "fileName")
          colnames(rg_tmp)[2:3] <- c("smin_new","smax_new")
        }else{
          rg_tmp <-
            signal_range[[x]] %>%
            data.frame(check.names = FALSE) %>%
            tibble::rownames_to_column(var = "fileName") %>%
            dplyr::mutate(smin_new = 0,.before = ".")
          colnames(rg_tmp)[3] <- "smax_new"
        }

        # rg_tmp <- signal_range[[x]] %>%
        #   data.frame() %>%
        #   tibble::rownames_to_column(var = "fileName")
        # colnames(rg_tmp)[2] <- "smax_new"

        # merge
        tmp <- rg_info %>% dplyr::filter(gene == names(signal_range)[x]) %>%
          left_join(.,rg_tmp,by = "fileName")
        # mutate(smax = smax_new)
      }) -> new_range

      # recover old range if is smax_new NA
      new_range <- new_range %>%
        mutate(smin_new = ifelse(is.na(smin_new),smin,smin_new),
               smax_new = ifelse(is.na(smax_new),smax,smax_new))
      # x = 1
      plyr::ldply(1:nrow(rg_info),function(x){
        tmp = rg_info[x,]
        tmp1 <- new_range %>% dplyr::filter(fileName == tmp$fileName & gene == tmp$gene)

        tmp <- tmp %>%
          mutate(smax = ifelse(nrow(tmp1) == 0,smax,tmp1$smax_new),
                 smin = ifelse(nrow(tmp1) == 0,smin,tmp1$smin_new))
      }) -> new_range

    }else{
      if(fixed_column_range == TRUE){
        plyr::ldply(1:length(unique(rg_info$gene)),function(x){
          tmp <- rg_info %>% dplyr::filter(gene == unique(rg_info$gene)[x]) %>%
            mutate(smax = max(smax),smin = min(smin))
        }) -> new_range
      }else{
        new_range <- rg_info
      }
    }

    # signal_layer_loop_params = list()
    # signal_layer_loop_params = list(max.height = 12)

    # add label yPos
    new_range <- new_range %>% mutate(rg_ypos = (smax - smin)*signal_range_pos[2] + smin)

    # re-assign loop height
    if(!is.null(signal_layer_loop_params$max.height)){
      plyr::ldply(1:nrow(new_range),function(x){
        tmp <- new_range[x,]
        if(tmp$track_type == "loop"){
          tmp$smax <- signal_layer_loop_params$max.height
        }
        return(tmp)
      }) -> new_range
    }

    # add range text label column
    new_range <- new_range %>%
      mutate(smax_value = ceiling(smax),
             smin_value = floor(smin),
             smax_label = paste("[",as.integer(floor(smin)),"-",as.integer(ceiling(smax)),"]",sep = ""))

    # ============================
    # whether reverse y axis
    # reverse_y_vars = NULL
    # reverse_y_vars = c("M1-CTCF","C1-CTCF","M1-H3K27ac","C1-H3K27ac")
    new_range <- new_range %>%
      mutate(yr_type = ifelse(fileName %in% reverse_y_vars,"reverse","identity"))

    # order
    new_range$fileName <- factor(new_range$fileName,levels = levels(tmp2$fileName))

    # layer
    # signal_range_label_params = list(color = "black",size = 4)
    range_label_layer <-
      # geom_text(data = new_range,
      #           aes(x = rg_xpos,y = rg_ypos,label = smax_label),
      #           color = signal_range_label_params$color,
      #           size = signal_range_label_params$size)
      # zplyr::geom_abs_text(data = new_range,
      #                      aes(xpos = signal_range_pos[1],
      #                          ypos = signal_range_pos[2],
      #                          label = smax_label),
      #                      color = signal_range_label_params$color,
      #                      size = signal_range_label_params$size)
      do.call(zplyr::geom_abs_text,
              modifyList(list(data = new_range,
                              aes(xpos = signal_range_pos[1],
                                  ypos = signal_range_pos[2],
                                  label = smax_label),
                              color = "black",
                              size = 4),signal_range_label_params))

    # mutate panel order
    new_range <- new_range %>% arrange(sample_group2,sample_group,
                                       fileName,
                                       gene_group2,gene_group,
                                       gene) %>%
      mutate(panel_num = 1:nrow(new_range))

    # trans panel y range
    tarns_pos_y_range <- unique(arrow.df[,c("gene","y")])

    # panel_range_layer
    iter_loop <- 1:(nrow(new_range) + length(add_facet_name)*(length(unique(tmp2$gene))))
    # if(is.null(Input_bed)){
    #   iter_loop <- iter_loop
    # }else{
    #   iter_loop <- iter_loop[-(length(unique(region.df$fileName)) + 1)]
    # }
    # x = 2
    panel_range_layer <- lapply(iter_loop, function(x){
      if(x <= nrow(new_range)){
        tmp <- new_range %>% dplyr::filter(panel_num == x)
        if(tmp$yr_type == "reverse"){
          sy <- scale_y_continuous(limits = c(tmp$smax_value,tmp$smin_value),trans = "reverse")
        }else{
          sy <- scale_y_continuous(limits = c(tmp$smin_value,tmp$smax_value))
        }
      }else{
        if(!is.null(Input_bed) & x %in%
           c((nrow(new_range) + 1):(nrow(new_range) + length(unique(Input_bed$sampleName))))){
          sy <- scale_y_continuous(limits = c(0,length(unique(Input_bed$sampleName)) + 1))
        }else{
          if(collapse_trans == TRUE){
            if(!is.null(Input_gene)){
              sy <- scale_y_continuous(limits = c(0,1 + 1.5))
            }else{
              sy <- scale_y_continuous(limits = c(0,2.5))
            }
          }else{
            if(!is.null(Input_gene)){
              if(!is.null(Input_bed)){
                index <- x - nrow(new_range) - length(unique(Input_bed$sampleName))
              }else{
                index <- x - nrow(new_range)
              }
              sy <- scale_y_continuous(limits = c(0,tarns_pos_y_range[index,]$y + 0.5))
            }else{
              tmp_p <- final_arrow_data %>% dplyr::select(gene,y) %>%
                group_by(gene) %>% summarise(y = max(y))
              sy <- scale_y_continuous(limits = c(0,tmp_p$y[x - nrow(new_range)] + 1.5))

            }
          }
        }
      }
      return(sy)
    })

    # panel x axis limit range
    # xlimit_range = list(c(1,2),c(3,4))
    iter_loop_x <- 1:(length(unique(tmp2$fileName))*length(unique(tmp2$gene)))
    if(!is.null(xlimit_range)){
      # check length
      if(length(xlimit_range) == 2 & !is.list(xlimit_range)){
        panel_x_range_layer <- lapply(iter_loop_x, function(x){
          sx <- scale_x_continuous(limits = c(xlimit_range[1],xlimit_range[2]))
          return(sx)
        })
      }else{
        panel_pos = matrix(iter_loop_x,ncol = length(unique(tmp2$gene)),byrow = T)
        panel_x_range_layer <-
          lapply(iter_loop_x, function(x){
            col_pos <- which(panel_pos == x, arr.ind = TRUE)[2]
            sx <- scale_x_continuous(limits = c(xlimit_range[[col_pos]]))
            return(sx)
          })
      }
    }else{
      panel_x_range_layer <- NULL
    }
  }else{
    # ==========================================================================
    range_label_layer <- NULL

    # trans panel y range
    tarns_pos_y_range <- unique(arrow.df[,c("gene","y")])

    iter_loop <- 1:(length(add_facet_name)*(length(unique(transcript.df$gene))))
    panel_range_layer <- lapply(iter_loop, function(x){
      if(collapse_trans == TRUE){
        sy <- scale_y_continuous(limits = c(0,2.5))
      }else{
        sy <- scale_y_continuous(limits = c(0,tarns_pos_y_range[x,]$y + 0.5))
      }
      return(sy)
    })

    # panel x axis limit range
    # xlimit_range = list(c(1,2),c(3,4))
    iter_loop_x <- 1:length(unique(tmp2$fileName))*length(unique(tmp2$gene))
    if(!is.null(xlimit_range)){
      # check length
      if(length(xlimit_range) == 2 & !is.list(xlimit_range)){
        panel_x_range_layer <- lapply(iter_loop_x, function(x){
          sx <- scale_x_continuous(limits = c(xlimit_range[1],xlimit_range[2]))
          return(sx)
        })
      }else{
        panel_pos = matrix(iter_loop_x,ncol = length(unique(tmp2$gene)),byrow = T)
        panel_x_range_layer <-
          lapply(iter_loop_x, function(x){
            col_pos <- which(panel_pos == x, arr.ind = TRUE)[2]
            sx <- scale_x_continuous(limits = c(xlimit_range[[col_pos]]))
            return(sx)
          })
      }
    }else{
      panel_x_range_layer <- NULL
    }
  }
  # ==============================================================================
  # 7_highlight mark region
  # ==============================================================================
  # higlight_region <- list(Actb = list(start = c(142904000),end = c(142904500)),
  #                   Myc = list(start = c(61986000,61989000),end = c(61986500,61989300)))
  #
  # higlight_col <- list(Actb = c("grey50"),
  #                      Myc = c("orange","green"))

  # higlight_region = NULL
  # higlight_col = NULL
  if(is.null(higlight_region) | is.null(higlight_col)){
    higlight_region_layer <- NULL
  }else{
    # x = 2
    plyr::ldply(1:length(higlight_region),function(x){
      tmp <- higlight_region[[x]]

      facet_info <- tmp2 %>%
        dplyr::filter(!(fileName %in% c("chrom")) & gene %in% names(higlight_region)[x]) %>%
        dplyr::select(fileName,gene,gene_group,gene_group2,sample_group,sample_group2) %>%
        distinct()

      res <- data.frame(xmin = tmp[["start"]],
                        xmax = tmp[["end"]],
                        ymin = -Inf,
                        ymax = Inf,
                        gene = names(higlight_region)[x],
                        col = higlight_col[[x]]) %>%
        left_join(.,facet_info,by = "gene",multiple = "all")
    }) -> hl_region

    # order
    hl_region$gene <- factor(hl_region$gene,levels = levels(tmp2$gene))

    # higlight_region_layer
    higlight_region_layer <-
      lapply(1:length(higlight_region), function(x){
        tmp <- hl_region %>% dplyr::filter(gene %in% names(higlight_region)[x])
        col_n <- unique(tmp$col)
        lapply(1:length(col_n), function(x){
          tmp_col <- tmp %>% dplyr::filter(col == col_n[x])
          geom_rect(data = tmp_col,
                    aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax),
                    fill = col_n[x],color = NA,
                    alpha = higlight_col_alpha,
                    show.legend = FALSE)
        }) -> layer_list
        return(layer_list)
      })
  }

  # ==============================================================================
  # 8_background color for each panel
  # ==============================================================================
  # background_color_region <- list(Actb = c("Input" = "grey90","IP-1"="yellow","IP-2"="purple"),
  #                                 Myc = c("Input" = "red","IP-1"="blue","IP-2"="green"))
  #
  # background_color_region <- list("chr5:65438603-65834288" = c("Input" = "grey90","IP-1"="yellow","IP-2"="purple"),
  #                                 "chr15:9263295-9785835" = c("Input" = "red","IP-1"="blue","IP-2"="green"))

  # background_color_region = NULL
  if(!is.null(background_color_region)){
    # x = 1
    plyr::ldply(1:length(background_color_region),function(x){
      tmp <- background_color_region[[x]]

      back_col_info <- tmp2 %>%
        dplyr::filter(gene %in% names(background_color_region)[x]) %>%
        dplyr::select(fileName,gene,gene_group,gene_group2,sample_group,sample_group2) %>%
        distinct() %>%
        mutate(xmin = -Inf,xmax = Inf,
               ymin = -Inf,ymax = Inf,
               col = tmp[fileName]) %>%
        dplyr::filter(col != "NA")
      return(back_col_info)

    }) -> background_region

    # background_region_layer
    background_region_layer <-
      lapply(1:nrow(background_region), function(x){
        tmp <- background_region[x,]
        geom_rect(data = tmp,
                  aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax),
                  fill = tmp$col,color = NA,
                  alpha = background_region_alpha,
                  show.legend = FALSE) -> layer_list
        return(layer_list)
      })
  }else{
    background_region_layer <- NULL
  }

  # ==============================================================================
  # 9_facet strips and labels settings
  # ==============================================================================
  # column_strip_setting_list = list(col = NA)
  # row_strip_setting_list = NULL
  # column_strip_text_setting_list = NULL
  # row_strip_text_setting_list = list(col = NA)
  # add_gene_region_label = FALSE

  # by_layer_x = FALSE
  # by_layer_y = FALSE

  # strip setting
  facet_strips <- ggh4x::strip_themed(
    background_x = do.call(ggh4x::elem_list_rect,
                           modifyList(list(colour = rep("white",2*length(unique(region.df$gene)))),
                                      column_strip_setting_list)),
    by_layer_x = by_layer_x,
    text_x = do.call(ggh4x::elem_list_text,
                     modifyList(list(),
                                column_strip_text_setting_list)),


    background_y = do.call(ggh4x::elem_list_rect,
                           modifyList(list(colour = rep("white",2*length(unique(region.df$gene)))),
                                      row_strip_setting_list)),
    by_layer_y = by_layer_y,
    text_y = do.call(ggh4x::elem_list_text,
                     modifyList(list(),
                                row_strip_text_setting_list)),
  )

  # whether supply with parames
  if(length(column_strip_setting_list) == 0 &
     length(column_strip_text_setting_list) == 0 &
     length(row_strip_setting_list) == 0 &
     length(row_strip_text_setting_list) == 0){
    facet_strips <- ggh4x::strip_nested()
  }else{
    facet_strips <- facet_strips
  }

  # gene facet label
  # add_gene_region_label = FALSE
  if(add_gene_region_label == TRUE){
    gene_chrom_df <-
      region.df %>% group_by(gene,seqnames) %>%
      summarise(start = min(start),end = min(end)) %>%
      mutate(seqnames = ifelse(startsWith(as.character(seqnames),"chr"),
                               seqnames,paste("chr",seqnames,sep = ""))) %>%
      mutate(label = paste(gene,"\n",seqnames,":",start,"-",end,sep = ""))

    new_column_label <- c(gene_chrom_df$label)
    names(new_column_label) <- gene_chrom_df$gene
  }else{
    new_column_label = NULL
  }

  # ==============================================================================
  # 10_facet settings
  # ==============================================================================
  i = if(is.null(gene_group_info)){FALSE}else{TRUE}
  j = if(is.null(gene_group_info2)){FALSE}else{TRUE}
  k = if(is.null(sample_group_info)){FALSE}else{TRUE}
  l = if(is.null(sample_group_info2)){FALSE}else{TRUE}


  if (!i & !j & !k & !l) {
    facet_var <- fileName~gene
  } else if (i & !j & !k & !l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is NULL, sample_group_info is NULL, sample_group_info2 is NULL\n")
    facet_var <- fileName~gene_group + gene
  } else if (!i & j & !k & !l) {
    # cat("gene_group_info is NULL, gene_group_info2 is not NULL, sample_group_info is NULL, sample_group_info2 is NULL\n")
    facet_var <- fileName~gene_group2 + gene
  } else if (!i & !j & k & !l) {
    # cat("gene_group_info is NULL, gene_group_info2 is NULL, sample_group_info is not NULL, sample_group_info2 is NULL\n")
    facet_var <- sample_group + fileName~gene
  } else if (!i & !j & !k & l) {
    # cat("gene_group_info is NULL, gene_group_info2 is NULL, sample_group_info is NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + fileName~gene
  } else if (i & j & !k & !l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is not NULL, sample_group_info is NULL, sample_group_info2 is NULL\n")
    facet_var <- fileName~gene_group2 + gene_group + gene
  } else if (i & !j & k & !l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is NULL, sample_group_info is not NULL, sample_group_info2 is NULL\n")
    facet_var <- sample_group + fileName~gene_group + gene
  } else if (i & !j & !k & l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is NULL, sample_group_info is NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + fileName~gene_group + gene
  } else if (!i & j & k & !l) {
    # cat("gene_group_info is NULL, gene_group_info2 is not NULL, sample_group_info is not NULL, sample_group_info2 is NULL\n")
    facet_var <- sample_group + fileName~gene_group2 + gene
  } else if (!i & j & !k & l) {
    # cat("gene_group_info is NULL, gene_group_info2 is not NULL, sample_group_info is NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + fileName~gene_group2 + gene
  } else if (!i & !j & k & l) {
    # cat("gene_group_info is NULL, gene_group_info2 is NULL, sample_group_info is not NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + sample_group + fileName~gene
  } else if (i & j & k & !l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is not NULL, sample_group_info is not NULL, sample_group_info2 is NULL\n")
    facet_var <- sample_group + fileName~gene_group2 + gene_group + gene
  } else if (i & j & !k & l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is not NULL, sample_group_info is NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + fileName~gene_group2 + gene_group + gene
  } else if (i & !j & k & l) {
    # cat("gene_group_info is not NULL, gene_group_info2 is NULL, sample_group_info is not NULL, sample_group_info2 is not NULL\n")
    facet_var <- sample_group2 + sample_group + fileName~gene_group + gene
  } else if (!i & j & k & l) {
    # cat("gene_group_info is NULL, gene_group")
    facet_var <- sample_group2 + sample_group + fileName~gene_group2 + gene
  }else{
    facet_var <- sample_group2 + sample_group + fileName~gene_group2 + gene_group + gene
  }


  facet_layer <- do.call(ggh4x::facet_nested,
                         modifyList(list(facet_var,
                                         scales = "free",
                                         nest_line = element_line(linetype = "solid"),
                                         independent = "y",
                                         switch = "y",
                                         strip = facet_strips,
                                         labeller = labeller(gene = new_column_label)),
                                    list()))

  # width and height of each panel
  # panel_size_setting = list()
  if(length(panel_size_setting) == 0){
    panel_size_layer <- NULL
  }else{
    panel_size_layer <- do.call(ggh4x::force_panelsizes,modifyList(
      list(respect = TRUE),panel_size_setting))
  }
  # ==============================================================================
  # 11_peaks facet plot
  # ==============================================================================
  # peaks_layer_params = list()
  if(!is.null(Input_bed)){
    peaks_layer <- do.call(geom_rect,
                           modifyList(list(data = bed.df.new,
                                           aes(xmin = start,xmax = end,
                                               ymin = ymin,ymax = ymax,
                                               fill = sampleName),
                                           color = NA),
                                      peaks_layer_params))
  }else{
    peaks_layer <- NULL
  }

  # ==============================================================================
  # 12_signal layer
  # ==============================================================================
  # signal_layer_bw_params = list()
  # signal_layer_loop_params = list()

  # bigwig layer
  # bw_geom = "rect"
  if(!is.null(Input_bw)){
    bw_data <- tmp2[which(tmp2$track_type == "bigwig"),]
    bw_data$start <- c(bw_data$start[1],bw_data$start[2:nrow(bw_data)]-1)

    signal_layer_geom_rect <- do.call(geom_rect,
                                      modifyList(
                                        list(data = bw_data,
                                             mapping = aes(xmin = start,xmax = end,
                                                           ymin = 0, ymax = score,
                                                           fill = fileName),
                                             color = NA,size = 0,
                                             show.legend = FALSE),
                                        signal_layer_bw_params))
  }else{
    signal_layer_geom_rect <- NULL
  }

  # loop layer
  if(!is.null(Input_loop)){
    if(Loop_curve_geom == "geom_arch"){
      signal_layer_geom_arch <- do.call(ggbio::geom_arch,
                                        modifyList(
                                          list(data = tmp2[which(tmp2$track_type == "loop"),],
                                               mapping = aes(x = start,xend = end,
                                                             height = score,
                                                             color = score),
                                               linewidth = 0.5,
                                               guide = guide_legend(FALSE)),
                                          signal_layer_loop_params))
    }else if(Loop_curve_geom == "geom_arch2"){
      signal_layer_geom_arch <- do.call(jjPlot::geom_arch2,
                                        modifyList(
                                          list(data = tmp2[which(tmp2$track_type == "loop"),],
                                               mapping = aes(x = start,xend = end,
                                                             y = 0,yend = score,
                                                             color = score)),
                                          signal_layer_loop_params))
    }else{
      message("Please supply 'geom_arch' or 'geom_arch2'!")
    }

  }else{
    signal_layer_geom_arch <- NULL
  }

  # heatmap layer
  if(!is.null(Input_hic)){
    tmp_ht_data <- tmp2[which(tmp2$track_type == "heatmap"),] %>%
      arrange(id)

    # signal_layer_heatmap_params = list()
    signal_layer_geom_polygon <- do.call(geom_polygon,
                                         modifyList(
                                           list(data = tmp_ht_data,
                                                mapping = aes(x = start,y = end,
                                                              group = id,fill = score),
                                                color = NA),
                                           signal_layer_heatmap_params))
  }else{
    signal_layer_geom_polygon <- NULL
  }

  # ============================================================================
  # junction data process
  if(!is.null(Input_junction)){
    # junc_layer_combined = FALSE
    if(junc_layer_combined == TRUE){
      range_bw <- new_range[which(new_range$track_type == "bigwig"),]
      junc_data <- tmp2[which(tmp2$track_type == "junction"),] %>%
        mutate(dist = end - start)

      # g = 1;f = 1
      file_name <- unique(junc_data$fileName)
      gene <- unique(junc_data$gene)
      plyr::ldply(seq_along(gene),function(g){
        tmp_g <- range_bw[which(range_bw$gene %in% gene[g]),]
        plyr::ldply(seq_along(file_name),function(f){
          match_bw_rg <- tmp_g[which(tmp_g$fileName %in% file_name[f]),]$smax_value
          tmp <- junc_data[which(junc_data$fileName %in% file_name[f] &
                                   junc_data$gene %in% gene[g]),]
          tmp$dist <- scales::rescale(tmp$dist,to = c(0.5*match_bw_rg,0.9*match_bw_rg))

          return(tmp)
        }) -> tmp1
        return(tmp1)
      }) -> junc_data
    }else{
      junc_data <- tmp2[which(tmp2$track_type == "junction"),] %>%
        mutate(dist = score)
    }

    # whether add junction band_line
    if(add_band_line == TRUE){
      total_count <- junc_data %>% group_by(gene,fileName) %>%
        summarise(sum_score = sum(score))

      # x = 1
      plyr::ldply(1:nrow(junc_data),function(x){
        tmp <- junc_data[x,]
        x_pos <- c(tmp$start,(tmp$start + tmp$end)/2,tmp$end)
        y_pos1 <- c(0,tmp$dist,0)

        # x_pos <- seq(tmp$start,tmp$end,length.out = 100)
        # y_pos1 <- c(0,tmp$dist,0)

        tmp_total <- total_count[which(total_count$gene == tmp$gene &
                                         total_count$fileName == tmp$fileName),]$sum_score

        # band_width = 0.5
        band_width_ratio <- band_width*tmp$dist*(tmp$score/tmp_total)

        y_pos2 <- c(0,tmp$dist - band_width_ratio,0)

        xy_coord1 <- data.frame(DescTools::DrawBezier(x = x_pos, y = y_pos1,
                                                      nv = 1000,plot = FALSE))
        xy_coord2 <- data.frame(DescTools::DrawBezier(x = x_pos, y = y_pos2,
                                                      nv = 1000,plot = FALSE))
        xy_coord2$x <- rev(xy_coord2$x)

        # res
        curve_band <- plyr::rbind.fill(xy_coord1,xy_coord2) %>%
          dplyr::mutate(seqnames = tmp$seqnames,
                        start = tmp$start,
                        end = tmp$end,
                        score = tmp$score,
                        fileName = factor(tmp$fileName,levels = levels(junc_data$fileName)),
                        track_type = factor(tmp$track_type,levels = levels(junc_data$track_type)),
                        gene = factor(tmp$gene,levels = levels(junc_data$gene)),
                        gene_group = factor(tmp$gene_group,levels = levels(junc_data$gene_group)),
                        gene_group2 = factor(tmp$gene_group2,levels = levels(junc_data$gene_group2)),
                        sample_group = factor(tmp$sample_group,levels = levels(junc_data$sample_group)),
                        sample_group2 = factor(tmp$sample_group2,levels = levels(junc_data$sample_group2)))

        curve_band$id = rep(x,2000)

        return(curve_band)
      }) -> band_junction_df

      # add geom_pologon
      signal_layer_junction <- do.call(geom_polygon,
                                       modifyList(
                                         list(data = band_junction_df,
                                              mapping = aes(x = x,y = y,
                                                            fill = fileName,
                                                            color = fileName,
                                                            group = id),
                                              color = NA,
                                              show.legend = FALSE),
                                         signal_layer_junction_params))

      # label data
      label_df <- band_junction_df %>%
        group_by(seqnames,start,end,fileName,track_type,gene,id,score) %>%
        summarise(y = max(y))

      signal_layer_junction_label <- do.call(geom_label,
                                             modifyList(
                                               list(data = label_df,
                                                    mapping = aes(x = (start + end)/2,
                                                                  y = y,
                                                                  label = score),
                                                    label.size = NA),
                                               signal_layer_junction_label_params))

    }else{
      # signal_layer_junction_params = list()
      # signal_layer_junction_label_params = list()
      if(Loop_curve_geom == "geom_arch"){
        signal_layer_junction <- do.call(ggbio::geom_arch,
                                         modifyList(
                                           list(data = junc_data,
                                                mapping = aes(x = start,xend = end,
                                                              color = fileName,
                                                              height = dist),
                                                linewidth = 0.5,
                                                show.legend = FALSE,
                                                guide = guide_legend(FALSE)),
                                           signal_layer_junction_params))
      }else if(Loop_curve_geom == "geom_arch2"){
        signal_layer_junction <- do.call(jjPlot::geom_arch2,
                                         modifyList(
                                           list(data = junc_data,
                                                mapping = aes(x = start,xend = end,
                                                              y = 0,yend = dist,
                                                              color = fileName),
                                                linewidth = 0.5,
                                                show.legend = FALSE),
                                           signal_layer_junction_params))
      }else{
        message("Please supply 'geom_arch' or 'geom_arch2'!")
      }


      signal_layer_junction_label <- do.call(geom_label,
                                             modifyList(
                                               list(data = junc_data,
                                                    mapping = aes(x = (start + end)/2,y = dist,
                                                                  label = score),
                                                    label.size = NA),
                                               signal_layer_junction_label_params))
    }

    #
    # ggplot(band_junction_df,aes(x,y,group = id)) +
    #   # geom_point() +
    #   geom_polygon(fill = "orange",color = NA) +
    #   geom_label(data = junc_data,
    #              mapping = aes(x = (start + end)/2,y = dist,
    #                            label = score),
    #              label.size = NA)

    # ggplot() +
    #   ggbio::geom_arch(data = junc_data,
    #                    aes(x = start,xend = end,
    #                        color = fileName,
    #                        height = dist)) +
    #   geom_label(data = junc_data,
    #              aes(x = (start + end)/2,y = dist,
    #                  label = score),
    #              label.size = NA) +
    # ggrepel::geom_label_repel(data = junc_data,
    #                           aes(x = (start + end)/2,y = dist,
    #                               label = score),
    #                           label.size = NA,
    #                           max.overlaps = 2000) +
    # theme_classic() +
    # facet_wrap(~fileName,scales = "free_y")
  }else{
    signal_layer_junction <- NULL
    signal_layer_junction_label <- NULL
  }


  # for test
  # ggplot() +
  #   geom_polygon(data = tmp_ht_data,
  #                mapping = aes(x = start,y = end,
  #                              group = id,fill = score),
  #                color = NA)
  #   scale_fill_gradient(low = 'grey90',high = 'red') +
  #   facet_wrap(~fileName,scales = "free")
  #   coord_equal()

  # ==============================================================================
  # 13_main plot
  # ==============================================================================
  # color settings
  if(is.null(sample_fill_col)){
    # sample_fill_col = rep("grey50",length(unique(tmp2$fileName)) - 1)
    sample_fill_col = suppressMessages(jjAnno::useMyCol(platte = "stallion",
                                                        n = length(unique(tmp2$fileName)) - 1))
  }else{
    sample_fill_col = sample_fill_col
  }

  # peaks color
  if(is.null(peak_fill_col)){
    peak_fill_col = suppressMessages(jjAnno::useMyCol(platte = "stallion",
                                                      n = length(unique(tmp2$fileName)) - 1))
  }else{
    peak_fill_col = peak_fill_col
  }

  # loop color
  if(is.null(loops_col)){
    loops_color_col = c("#5BC0F8","#FF597B")
  }else{
    loops_color_col = loops_col
  }

  # heatmap color
  if(is.null(heatmap_fill_col)){
    heatmap_fill_col = RColorBrewer::brewer.pal(n = 9,name = "Greens")
  }else{
    heatmap_fill_col = heatmap_fill_col
  }

  # exon colors
  if(is.null(trans_fill_col)){
    trans_fill_col = RColorBrewer::brewer.pal(n = 9,name = "Paired")
  }else{
    trans_fill_col = trans_fill_col
  }

  if(draw_chromosome == TRUE){
    useless_col = rep("white",2)
  }else{
    useless_col = rep("white",1)
  }

  # ============================================================
  # whether show y axis ticks
  if(show_y_ticks == TRUE){
    range_label_layer <- NULL
    axis.text.y <- element_text()
    axis.ticks.y <- element_line()
  }else{
    axis.text.y <- element_blank()
    axis.ticks.y <- element_blank()
  }

  # ============================================================
  # combine all layers
  pmain <-
    ggplot() +
    background_region_layer +
    # loop layer
    signal_layer_geom_arch +
    scale_color_gradientn(colors = loops_color_col,
                          name = "loop_distance\n(10kb)") +
    # bigwig layer
    ggnewscale::new_scale_fill() +
    # ggnewscale::new_scale_color() +
    signal_layer_geom_rect +
    scale_fill_manual(values = c(sample_fill_col,useless_col),
                      name = "") +
    # scale_color_manual(values = ggplot2::alpha(c(sample_fill_col,useless_col),
    #                                            alpha = ifelse(!is.null(signal_layer_bw_params$alpha),
    #                                                           signal_layer_bw_params$alpha,1)),
    #                    name = "") +
    # heatmap layer
    ggnewscale::new_scale_fill() +
    signal_layer_geom_polygon +
    scale_fill_gradientn(colors = heatmap_fill_col,name = "contact\nfrequency") +
    # bed layer
    ggnewscale::new_scale_fill() +
    peaks_layer +
    scale_fill_manual(values = c(peak_fill_col,useless_col),name = "") +
    # junction layer
    ggnewscale::new_scale_colour() +
    signal_layer_junction +
    scale_color_manual(values = c(sample_fill_col,useless_col),
                       name = "") +
    signal_layer_junction_label +
    # other layers
    text_layer +
    seg_arrow_layer +
    # transcript layer
    trans_arrow_layer +
    ggnewscale::new_scale_fill() +
    trans_struct_layer +
    scale_fill_manual(values = trans_fill_col,name = "exon type") +
    # gene label layer
    gene_label_layer +
    higlight_region_layer +
    theme_bw(base_size = base_size) +
    facet_layer +
    ggh4x::facetted_pos_scales(y = panel_range_layer,
                               x = panel_x_range_layer) +
    range_label_layer +
    panel_size_layer +
    theme(panel.grid = element_blank(),
          panel.spacing.y = unit(panel.spacing[2],"lines"),
          panel.spacing.x = unit(panel.spacing[1],"lines"),
          strip.text.y.left = element_text(angle = 0,face = "bold",hjust = 1),
          strip.text.x = element_text(face = "bold.italic"),
          axis.text.y = axis.text.y,
          axis.ticks.y = axis.ticks.y,
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = "black")) +
    xlab("") + ylab("")

  # transcript.df %>% group_by(gene) %>%
  #   summarise(start = min(start),end = max(end))
  # ==============================================================================
  # 14_draw chromosome panel
  # ==============================================================================
  if(draw_chromosome == TRUE){
    annoChr <- lapply(1:length(unique(segment.df$gene)), function(x){
      tmp <- segment.df %>% dplyr::filter(gene == unique(segment.df$gene)[x]) %>%
        mutate(fileName = "chrom",
               seqnames = if_else(startsWith(as.character(seqnames),"chr"),
                                  seqnames,paste("chr",seqnames,sep = "")))

      # order
      t_levels <- if(!is.null(Input_gene)){
        levels(segment.df$fileName)
      }else{
        levels(tmp2$fileName)
      }
      tmp$fileName <- factor(tmp$fileName,levels = t_levels)

      # chromosome grob
      # draw_chromosome_params = list(karyotype_file = NULL,centromereTri_file = NULL,density_file = NULL)
      if(!is.null(Input_gene)){
        chromosomes = tmp$seqnames
        zoom_region = c(tmp$start,tmp$end)
      }else{
        region_tmp <- region.df[which(region.df$gene %in% unique(segment.df$gene)[x]),]
        chromosomes = unique(tmp$seqnames)
        zoom_region = c(min(region_tmp$start),max(region_tmp$start))
      }

      # check whether set xlimit
      if(!is.null(xlimit_range)){
        if(!is.list(xlimit_range)){
          zoom_region = xlimit_range
        }else{
          zoom_region = xlimit_range[[x]]
        }
      }

      pc <-
        # drawChromosome(karyotype_file = karyotype.df,
        #                centromereTri_file = centromere.triangle,
        #                density_file = density.df,
        #                chromosomes = tmp$seqnames,
        #                zoom_region = c(tmp$start,tmp$end),
        #                facet_params = list(ncol = 1,strip_position = "bottom"),
        #                add_regionLen = FALSE) +
        do.call(drawChromosome,modifyList(list(chromosomes = chromosomes,
                                               zoom_region = zoom_region,
                                               facet_params = list(ncol = 1,strip.position = "bottom"),
                                               add_regionLen = FALSE),
                                          draw_chromosome_params)) +
        scale_x_continuous(expand = c(0,0)) +
        theme_void() +
        theme(strip.text = element_blank(),
              plot.margin = margin(0,0,0,0,"cm"))

      # inset to panels
      annotation_custom2(grob = ggplotGrob(pc),
                         data = tmp,
                         xmin = zoom_region[1], xmax = zoom_region[2],
                         ymin = -Inf, ymax = Inf)
    })

    # final plot
    pfinal <- pmain + annoChr
  }else{
    pfinal <- pmain
  }

  # ==============================================================================
  # 15_whether remove panel borders
  # ==============================================================================
  # remove_chrom_panel_border = FALSE
  # remove_all_panel_border = FALSE
  if(remove_chrom_panel_border == TRUE | remove_all_panel_border == TRUE){
    g <- ggplotGrob(pfinal)

    # i = 10
    # remove_all_panel_border = TRUE
    if(is.null(Input_gene)){
      col_num <- length(query_region[[1]])
    }else{
      col_num <- length(Input_gene)
    }

    # defnie panels
    if(remove_all_panel_border == TRUE){
      panel_num <- 2:(col_num*(length(unique(region.df$fileName)) + 2) + 1)
    }else if(remove_chrom_panel_border == TRUE){
      end = col_num*(length(unique(region.df$fileName)) + 2) + 1
      start = end - col_num + 1
      panel_num <- start:end
    }else{
      message("Should not be both TRUE!")
    }

    # grobs editting
    for (i in panel_num) {
      grobs_border <- grid::grid.ls(g$grobs[[i]],print = FALSE)
      panel_boder_name <- grobs_border[["name"]][length(grobs_border[["name"]])]
      g$grobs[[i]]$children[[panel_boder_name]]$gp$col <- NA
    }
    grid::grid.newpage()
    grid::grid.draw(g)
  }else{
    pfinal
  }

}
