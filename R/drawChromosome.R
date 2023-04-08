#' Draw Chromosome Karyotype Plot
#'
#' @author JunZhang
#'
#' This function creates a chromosome karyotype plot from the given input files
#' and parameters.
#'
#' @param ideogram_obj The output for "catchIdeoData" function.
#' @param karyotype_file Path to the karyotype file. The file is expected to be
#' in tab-delimited format with columns chr, start, end, and name representing
#' chromosome name, start and end positions, and chromosome band name respectively.
#' @param centromereTri_file Path to the centromere triangle file. The file is
#' expected to be in tab-delimited format with columns chr, start, and end
#' representing chromosome name, start and end positions of the centromere
#' triangle respectively.
#' @param density_file Path to the density file. The file is expected to be in
#' tab-delimited format with columns chr, pos, and density representing chromosome
#' name, position and density value.
#' @param chromosomes A character vector indicating the chromosomes to be
#' displayed in the plot.
#' @param orders A named integer vector indicating the order of chromosomes to be
#' displayed in the plot.
#' @param density_col A character vector of length 2 indicating the colors to be
#' used for the density plot.
#' @param centromere_col A character value indicating the color to be used for
#' the centromere triangles.
#' @param zoom_region A numeric vector of length 2 indicating the region to be
#' zoomed in. Default is NULL (no zoom).
#' @param zoom_pos A character value indicating the location of the zoomed region.
#' Possible values are "top" or "bottom".
#' @param zoom_col A character vector of length 2 indicating the colors to be
#' used for the zoomed region.
#' @param zoom_relheight A numeric value indicating the relative height of the
#' zoomed region with respect to the main plot.
#' @param facet_params List of parameters for facet_wrap function.
#' @param add_regionLen Whether add region length information.
#' @param ... Additional parameters to be passed to the theme function.
#'
#' @importFrom ggplot2 ggplot geom_segment geom_rect geom_polygon aes scale_fill_manual
#' theme_bw theme xlab ylab element_blank element_text scale_fill_gradient arrow unit
#' geom_text
#' @importFrom dplyr filter group_by summarize select
#' @importFrom plyr ldply
#' @export
drawChromosome <- function(ideogram_obj = NULL,
                           karyotype_file = NULL,
                           centromereTri_file = NULL,
                           density_file = NULL,
                           chromosomes = NULL,
                           orders = NULL,
                           density_col = c("grey90","black"),
                           centromere_col = "#F00000",
                           zoom_region = NULL,
                           zoom_pos = "top",
                           zoom_col = c("grey95","grey60"),
                           zoom_relheight = 1,
                           facet_params = list(),
                           add_regionLen = TRUE,
                           ...){
  # ============================================================================
  # 1_basic plot
  # ============================================================================
  if(is.null(ideogram_obj)){
    # order
    if(!is.null(orders)){
      karyotype_file$chr <- factor(karyotype_file$chr,levels = orders)
      centromereTri_file$chr <- factor(centromereTri_file$chr,levels = orders)
      density_file$chr <- factor(density_file$chr,levels = orders)
    }

    # filter
    if(!is.null(chromosomes)){
      karyotype_file <- karyotype_file %>% dplyr::filter(chr %in% chromosomes)
      centromereTri_file <- centromereTri_file %>% dplyr::filter(chr %in% chromosomes)
      density_file <- density_file %>% dplyr::filter(chr %in% chromosomes)
    }

    # plot
    pmain <-
      ggplot() +
      geom_rect(data = density_file,
                aes(xmin = start,xmax = end,
                    ymin = 0,ymax = 1,
                    fill = density),
                color = NA,show.legend = FALSE) +
      geom_rect(data = karyotype_file,
                aes(ymin = 0,ymax = 1,
                    xmin = start,xmax = end),
                color = "black",fill = "transparent") +
      # triangle graph
      geom_polygon(data = centromereTri_file,
                   aes(x = x,y = y),
                   color = NA,fill = centromere_col) +
      scale_fill_gradient(low = density_col[1],high = density_col[2])

  }else{
    # filter target chromosomes
    if(!is.null(chromosomes)){
      plot_df <- ideogram_obj$plot_df %>% dplyr::filter(chr %in% chromosomes)
      border_df <- ideogram_obj$border_df %>% dplyr::filter(chr %in% chromosomes)

      if(!is.null(ideogram_obj$acen_plot_df)){
        acen_plot_df <- ideogram_obj$acen_plot_df %>% dplyr::filter(chr %in% chromosomes)
      }else{
        acen_plot_df <- NULL
      }
    }else{
      plot_df <- ideogram_obj$plot_df
      border_df <- ideogram_obj$border_df
      acen_plot_df <- ideogram_obj$acen_plot_df
    }

    # check data
    if(!is.null(ideogram_obj$acen_plot_df)){
      centromere_layer <- geom_polygon(data = acen_plot_df,
                                       aes(x = xpos,y = ypos,group = type),
                                       color = "black",fill = "darkred")
    }else{
      centromere_layer <- NULL
    }

    # colors
    col_df <- plot_df %>% dplyr::select(gieStain,gieStainCol) %>%
      unique()
    col_c <- col_df$gieStainCol
    names(col_c) <- col_df$gieStain

    # plot
    pmain <-
      ggplot() +
      geom_rect(data = plot_df,
                aes(xmin = chromStart,xmax = chromEnd,
                    ymin = 0,ymax = 1,
                    fill = gieStain),show.legend = FALSE) +
      geom_rect(data = border_df,
                aes(xmin = xmin,xmax = xmax,
                    ymin = 0,ymax = 1),
                fill = NA,color = "black") +
      centromere_layer +
      scale_fill_manual(values = col_c)
  }

  # plot layout
  playout <- pmain +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold"),
          ...
    ) +
    ylab("") + xlab("") +
    do.call(ggplot2::facet_wrap,modifyList(list(~chr,
                                                ncol = 1,
                                                scales = "free_x",
                                                strip.position = "left"),
                                           facet_params))
  # ============================================================================
  # 2_whether add zoom track
  # ============================================================================
  if(is.null(ideogram_obj)){
    zoom_input <- karyotype_file
  }else{
    zoom_input <- plot_df %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(start = 0,end = max(chromEnd))
  }

  # loop
  # x = 1
  if(!is.null(zoom_region)){
    plyr::ldply(1:nrow(zoom_input),function(x){
      tmp <- zoom_input[x,]

      # if supply multiple positions
      if(!is.list(zoom_region)){
        xmin <- zoom_region[1]
        xmax <- zoom_region[2]
      }else{
        xmin <- zoom_region[[tmp$chr]][1]
        xmax <- zoom_region[[tmp$chr]][2]
      }

      # zoom position
      if(zoom_pos != "top"){
        xPos = c(tmp$start,xmin,xmax,tmp$end)
        yPos = c(-1*zoom_relheight,0,0,-1*zoom_relheight)
      }else{
        xPos = c(xmin,tmp$start,tmp$end,xmax)
        yPos = c(1,1+zoom_relheight,1+zoom_relheight,1)
      }

      # get data
      trapezoid.df <- createTrapezoid(xPos = xPos,yPos = yPos) %>%
        dplyr::mutate(chr = tmp$chr)

      return(trapezoid.df)
    }) -> target.zoom


    # higlight region
    if(!is.list(zoom_region)){
      hdf <- data.frame(chr = zoom_input$chr,
                        start = zoom_region[1],
                        end = zoom_region[2])
    }else{
      hdf <- t(data.frame(zoom_region)) %>% data.frame() %>%
        tibble::rownames_to_column(var = "chr")
      colnames(hdf)[2:3] <- c("start","end")
    }

    # sgement data
    segdf <- hdf %>% dplyr::left_join(.,zoom_input,by = "chr") %>%
      dplyr::mutate(label = paste(round(abs(end.x - start.x)/10^4,digits = 1),"kb",sep = " "),
                    segy = dplyr::if_else(zoom_pos == "top",
                                          1 + zoom_relheight + 0.1,-1*zoom_relheight - 0.1))


    # add zoom track
    pzoom <-
      playout +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = target.zoom,
                   aes(x = x, y = y,group = id,fill = id),
                   colour = NA,linewidth = 0,
                   show.legend = FALSE) +
      scale_fill_gradient(low = zoom_col[1],high = zoom_col[2]) +
      # add highlit region
      geom_rect(data = hdf,
                aes(xmin = start,xmax = end,
                    ymin = -0.05,ymax = 1.05),
                fill = "red",color = "red",size = 1)

    # add segment table
    if(add_regionLen == TRUE){
      pres <- pzoom +
        # add segment
        geom_segment(data = segdf,
                     aes(x = 0,xend = end.y,
                         y = segy,yend = segy),
                     color = "black",linewidth = 0.5,
                     arrow = arrow(type = "open",ends = "both",
                                   angle = 30,length = unit(2.5,"mm"))) +
        geom_text(data = segdf,aes(x = end.y/2,y = segy + 0.3,label = label))
    }else{
      pres <- pzoom
    }
  }else{
    pres <- playout
  }

  return(pres)
}


# =========================================================================================
# test data
# =========================================================================================

#' This is a test data for this package
#' test data describtion
#'
#' @name hg38_obj
#' @docType data
#' @author JunZhang
"hg38_obj"

#' This is a test data for this package
#' test data describtion
#'
#' @name hg19_obj
#' @docType data
#' @author JunZhang
"hg19_obj"

#' This is a test data for this package
#' test data describtion
#'
#' @name mm10_obj
#' @docType data
#' @author JunZhang
"mm10_obj"
