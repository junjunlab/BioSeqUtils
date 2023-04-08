globalVariables(c("browserSession", "chrom", "chromEnd", "chromStart", "genome" , "gieStain",
                  "idioLabelColors", "modifyList", "pos", "scale_fill_manual", "ucscSchema",
                  "ucscTableQuery", "xmax", "xmin", "xpos", "ypos","genome<-", "getTable",
                  "gieStainCol"))

#' Download and process cytoband ideogram data
#'
#' @author JunZhang
#'
#' @param genome The genome assembly version (default is "hg19")
#' @param chromosome The chromosome number to filter by (default is NULL)
#' @return A list containing the processed ideogram data
#'
#' @details This function downloads cytoband ideogram data from the UCSC genome browser
#' and assigns colors to different regions based on their densities. It then prepares the
#' data for plotting by creating border and centromere connection data, as well as a plot
#' dataframe. The resulting data is returned as a list.
#'
#' @description
#' * "gneg": 染色体中的负色区域。这些区域通常不含有基因，且在染色过程中呈现较浅的颜色.
#' * "gpos25": 染色体中呈现较深颜色的区域，具有25%的染色强度。这些区域通常含有基因.
#' * "gpos50": 染色体中呈现较深颜色的区域，具有50%的染色强度。这些区域通常含有基因.
#' * "gpos75": 染色体中呈现较深颜色的区域，具有75%的染色强度。这些区域通常含有基因.
#' * "gpos100": 染色体中呈现最深颜色的区域，具有100%的染色强度。这些区域通常含有基因.
#' * "acen": 染色体的着丝点。这是染色体中的一个特殊区域，负责染色体在细胞分裂过程中的移动与排列.
#' * "gvar": 染色体中的变异区域。这些区域在不同个体或种群之间具有较高的基因变异率.
#' * "stalk": 染色体的柄区。这是染色体中的一个狭窄区域，连接着着丝点和其他染色质区域.
#'
#' @examples
#' \dontrun{ideogram_data <- catchIdeoData(genome = "hg19")}
#'
#' @export
catchIdeoData <- function(genome = "hg19",
                          chromosome = NULL){
  options(warn = -1)
  options(dplyr.summarise.inform = FALSE)
  # ==========================================
  # 1_down load data
  # ==========================================
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, "cytoBandIdeo")
  ucscSchema(query)
  ideonew <- getTable(query)

  # ==========================================
  # 2_assign colors to different regions
  # ==========================================

  # setting band colors
  # gneg & gposN(where N is an integer [1,100] ) prepresent densities.
  # gneg has density 0
  if(startsWith(genome,"hg")){
    idioLabelColors <- c(gneg = "white", gpos25 = "gray75", gpos50 = "gray50",
                         gpos75 = "gray25", gpos100 = "black",acen = "darkred",
                         gvar = "black", stalk = "#19A7CE")
  }else{
    idioLabelColors <- c(gneg = "white", gpos33 = "gray67", gpos66 = "gray34",
                         gpos75 = "gray25", gpos100 = "black")
  }

  ideonew$gieStainCol <- idioLabelColors[match(ideonew$gieStain,names(idioLabelColors))]

  # filter chromosomes
  if(is.null(chromosome)){
    targetChrom <- paste("chr",c(1:22,"X","Y"),sep = "")
    ideonew <- ideonew[which(ideonew$chrom %in% targetChrom),]
  }else{
    ideonew <- ideonew[which(ideonew$chrom %in% targetChrom),]
  }

  gieStain_type <- unique(ideonew$gieStain)
  # ==========================================
  # 3_parepare data for plot
  # ==========================================
  acen_df <- ideonew[which(ideonew$gieStain == "acen"),] %>%
    group_by(chrom) %>%
    summarise(start = min(chromStart),end = max(chromEnd))

  # borders
  # x = 1
  if("acen" %in% gieStain_type){
    clean_df <- ideonew %>% filter(gieStain != "acen")
    plyr::ldply(1:nrow(clean_df),function(x){
      tmp <- clean_df[x,]
      acen_tmp <- acen_df[which(acen_df$chrom %in% tmp$chrom),]
      tmp <- tmp %>% mutate(pos = ifelse(chromEnd <= acen_tmp$start,"left","right"))
      return(tmp)
    }) -> border_df

    border_df <- border_df %>%
      group_by(chrom,pos) %>%
      summarise(xmin = min(chromStart),xmax = max(chromEnd))
  }else{
    border_df <- ideonew %>%
      group_by(chrom) %>%
      summarise(xmin = min(chromStart),xmax = max(chromEnd))
  }

  border_df <- border_df %>% rename(chr = chrom)

  # centromere connection data
  if("acen" %in% gieStain_type){
    # x = 1
    plyr::ldply(1:nrow(acen_df),function(x){
      tmp <- acen_df[x,]
      # get polygon position
      df <- createPairPolygon(start = tmp$start,end = tmp$end,open = 0.5)
      df$ypos <- df$ypos + 0.5
      df$chrom <- tmp$chrom
      return(df)
    }) -> acen_plot_df
    acen_plot_df <- acen_plot_df %>% rename(chr = chrom)
  }else{
    acen_plot_df <- NULL
  }


  # plot df
  plot_df <- ideonew %>% filter(gieStain != "acen") %>% rename(chr = chrom)

  # output
  ideogram_data = list(plot_df = plot_df,border_df = border_df,acen_plot_df = acen_plot_df)
  return(ideogram_data)
}
