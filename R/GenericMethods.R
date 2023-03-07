globalVariables(c("cdslen", "end", "start", 'tga', "tlen", ".",
                  "cdsed", "geneId", "geneName", "selecType",
                  "sep", "topN", "transId","cdsst","tname",
                  '5UTR', 'CDS', 'exon', 'gene_id', 'gene_name',
                  'gtype', 'transcript_id', 'type', 'typelen',
                  'width','mytest', 'n', 'seqnames', 'strand'))

# ==============================================================================
# setGeneric
# ==============================================================================
#' Define a generic method for showing objects
#'
#' This function defines a generic method for showing objects, which can be
#' implemented for different classes of objects.
#'
#' @title myShow
#' @param object An object to myShow.
#'
#' @return None.
#' @export
setGeneric("myShow",function(object) standardGeneric("myShow"))


#' filterID method for GenomeGTF objects
#'
#' @title filterID method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("filterID",function(object,...) standardGeneric("filterID"))


# setGeneric("filterRepTrans",function(object,...) standardGeneric("filterRepTrans"))


#' getLongName method for GenomeGTF objects
#'
#' @title getLongName method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("getTransInfo",function(object,...) standardGeneric("getTransInfo"))


#' getFeatureFromGenome method for GenomeGTF objects
#'
#' @title getFeatureFromGenome method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("getFeatureFromGenome",function(object,...) standardGeneric("getFeatureFromGenome"))


#' getIntronInfo method for GenomeGTF objects
#'
#' @title getIntronInfo method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("getIntronInfo",function(object,...) standardGeneric("getIntronInfo"))


# setGeneric("superExtract",function(object,...) standardGeneric("superExtract"))

# ==============================================================================
# setMethod
# ==============================================================================
setMethod("show",
          signature(object = "GenomeGTF"),
          function(object){
            cat("## GenomeGTF object for Extracting sequences.\n")
            if(is.null(object@gtf)){
              cat("## GTF file is NULL.\n")
            }else{
              cat("## GTF file is loaded.\n")
              cat(paste0("## GTF path: ",object@gtfPath,".\n",collapse = ""))
            }

            if(is.null(object@genome)){
              cat("## genome file is NULL.\n")
            }else{
              cat("## genome file is loaded.\n")
              cat(paste0("## Genome path: ",object@genomePath,".\n",collapse = ""))
            }

            if(is.null(object@representTrans)){
              cat("## representTrans file is NULL.\n")
            }else{
              cat("## representTrans file is loaded.\n")
            }

            # cat("## representTrans slot is NULL.\n")
            cat("## intron slot is NULL.\n")
          })

#' Show method for GenomeGTF object
#'
#' This function displays information about a GenomeGTF object, including the
#' status of its GTF and genome files.
#'
#' @param object An object of class GenomeGTF.
#'
#' @return Nothing is returned; this function prints information to the console.
#'
#' @examples
#' # Show information about the object
#' \dontrun{show(GenomeGTF)}
#'
#' @export
setMethod("myShow",
          signature(object = "GenomeGTF"),
          function(object){
            cat("## GenomeGTF object for Extracting sequences.\n")
            if(is.null(object@gtf)){
              cat("## GTF file is NULL.\n")
            }else{
              cat("## GTF file is loaded.\n")
            }

            if(is.null(object@genome)){
              cat("## genome file is NULL.\n")
            }else{
              cat("## genome file is loaded.\n")
            }

            if(is.null(object@representTrans)){
              cat("## representTrans file is NULL.\n")
            }else{
              cat("## representTrans file is loaded.\n")
            }

            # cat("## representTrans slot is NULL.\n")
            cat("## intron slot is NULL.\n")
          })

# filterID
# filter information in annotation file with given geneName,geneId or transId
# and return a filtered data frame.
#' Filter genome annotation information by gene name, gene ID, or transcript ID
#'
#' This function filters the genome annotation information in a GenomeGTF object
#' based on the provided gene name, gene ID, or transcript ID.
#'
#' @param object A GenomeGTF object.
#' @param geneName A character vector of gene names to filter by. Default is NULL.
#' @param geneId A character vector of gene IDs to filter by. Default is NULL.
#' @param transId A character vector of transcript IDs to filter by. Default is NULL.
#'
#' @return A data.frame containing the filtered genome annotation information
#' @method filterID GenomeGTF
#' @export
setMethod("filterID",
          signature(object = "GenomeGTF"),
          function(object,
                   geneName = NULL,geneId = NULL,transId = NULL){
            # load gtf
            gtf <- as.data.frame(object@gtf)

            if(!is.null(geneName) & is.null(geneId) & is.null(transId)){
              ginfo <- gtf[which(gtf$gene_name %in% geneName),]
            }else if(is.null(geneName) & !is.null(geneId) & is.null(transId)){
              ginfo <- gtf[which(gtf$gene_id %in% geneId),]
            }else if(is.null(geneName) & is.null(geneId) & !is.null(transId)){
              ginfo <- gtf[which(gtf$transcript_id %in% transId),]
            }else if(!is.null(geneName) & is.null(geneId) & !is.null(transId)){
              ginfo <- gtf[which(gtf$gene_name == geneName & gtf$transcript_id == transId),]
            }else if(is.null(geneName) & !is.null(geneId) & !is.null(transId)){
              ginfo <- gtf[which(gtf$gene_id == geneId & gtf$transcript_id == transId),]
            }else{
              message("Please choose again.")
            }

            # output
            return(ginfo)
          })



# setMethod("filterRepTrans",
#           signature(object = "GenomeGTF"),
#           function(object,
#                    geneName = NULL,
#                    geneId = NULL,
#                    selecType = c("lcds","lt"),
#                    # 0 for all selections
#                    topN = 1,sep = "|"){
#             # doing now
#             selecType <- match.arg(selecType)
#
#             # get data
#             if(as.character(class(object)) == "GenomeGTF"){
#               gtf <- as.data.frame(object@gtf)
#
#               # choose type
#               if(!is.null(geneName) & is.null(geneId)){
#                 tga <- gtf[which(gtf$gene_name %in% geneName),]
#               }else if(is.null(geneName) & !is.null(geneId)){
#                 tga <- gtf[which(gtf$gene_id %in% geneId),]
#               }else if(!is.null(geneName) & !is.null(geneId)){
#                 tga <- gtf[which(gtf$gene_name == geneName & gtf$gene_id == geneId),]
#               }
#
#             }else if(as.character(class(object)) == "data.frame"){
#               tga <- object[which(object$gene_id %in% geneId),]
#             }
#
#             gid = unique(tga$gene_id)
#
#             # progress bar
#             pb <- progress::progress_bar$new(
#               format = 'filterRepTrans is running [:bar] :percent in :elapsed',
#               total = length(gid), clear = FALSE, width = 80
#             )
#
#             # parallel calculation
#             # options(future.globals.maxSize = memorySize*1000*1024^2)
#             # future::plan("future::multisession")
#
#             # loop
#             plyr::ldply(1:length(gid),function(x){
#               # future.apply::future_lapply(1:length(gid),function(x){
#               pb$tick()
#
#               tg <- tga[which(tga$gene_id == gid[x]),]
#
#               # calculate transcript/CDS length
#               tida <- unique(tg$transcript_id)
#               plyr::ldply(1:length(tida),function(x){
#                 # future.apply::future_lapply(1:length(tida),function(x){
#
#                 tmp.info <- tg[which(tg$transcript_id == tida[x]),]
#                 t.type <- unique(tmp.info$type)
#
#                 # exon info
#                 tmp <- tmp.info[which(tmp.info$type == "exon"),]
#
#                 # check transcript whether has CDS
#                 if("CDS" %in% t.type){
#                   tmp2 <- tmp.info[which(tmp.info$type == "CDS"),]
#                   cdslen = sum(tmp2$width)
#                 }else{
#                   cdslen = 0
#                 }
#
#                 # return data.frame
#                 tranLength <- data.frame(gname = unique(tg$gene_name),
#                                          gid = unique(tg$gene_id),
#                                          tid = tida[x],
#                                          tlen = sum(tmp$width),
#                                          cdslen = cdslen)
#
#                 return(tranLength)
#               }) -> lenInfo
#               # }) %>% do.call("rbind",.) %>% data.frame() -> lenInfo
#
#               # choose type
#               if(selecType == "lt"){
#                 rankTran <- lenInfo %>%
#                   dplyr::arrange(dplyr::desc(tlen),dplyr::desc(cdslen))
#               }else if(selecType == "lcds"){
#                 rankTran <- lenInfo %>%
#                   dplyr::arrange(dplyr::desc(cdslen),dplyr::desc(tlen))
#               }else{
#                 message("Please choose 'lt' or 'lcds'!")
#               }
#
#               # get long name
#               rankTran.name <- getLongName(object = object,transId = rankTran$tid,sep = sep)
#               # rankTran$tname <- rankTran.name$tname
#
#               rankTran$tname <- ifelse(rankTran$cdslen == 0,
#                                        paste(rankTran.name$tname,"NC",sep = sep),
#                                        paste(rankTran.name$tname,"CD",sep = sep))
#
#               # return length info
#               if(topN == 0){
#                 return(rankTran)
#               }else{
#                 plyr::ldply(unique(rankTran$gid),function(x){
#                   rankTran <- rankTran[which(rankTran$gid == x),] %>%
#                     dplyr::slice_head(n = as.numeric(topN))
#                   return(rankTran)
#                 }) -> rankTran
#                 return(rankTran)
#               }
#
#               Sys.sleep(0.05)
#             }) -> rep.info
#             # }) %>% do.call("rbind",.) %>% data.frame() -> rep.info
#
#             return(rep.info)
#           })


# getTransInfo
# calculate CDS absolute start position and stop position for transcript.
#'Get long name of transcripts in GenomeGTF object
#'
#' This function extracts information about transcripts in a GenomeGTF object
#' and returns a data frame with their long names.
#'
#' @param object A GenomeGTF object.
#' @param geneName Character vector of gene names to filter by. Default is NULL.
#' @param geneId Character vector of gene IDs to filter by. Default is NULL.
#' @param transId Character vector of transcript IDs to filter by. Default is NULL.
#' @param selecType A character vector specifying the selection criteria.
#' Must be either "lt" for longest transcript or "lcds" for longest CDS. Default
#' is "lcds".
#' @param topN An integer indicating the number of longest transcripts to retain.
#' Default is 1. Set to 0 to retain all transcripts.
#' @param sep Separator character to use in the long name. Default is "|".
#'
#' @return A data frame.
#'
#' @importFrom dplyr mutate desc select arrange slice_head group_by
#' @method getTransInfo GenomeGTF
#' @export
setMethod("getTransInfo",
          signature(object = "GenomeGTF"),
          function(object,
                   geneName = NULL,geneId = NULL,transId = NULL,
                   selecType = c("lcds","lt"),
                   topN = 1,sep = "|"){
            # load GTF
            ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)

            # ===============================================================================
            # recode
            leninfo <- ginfo[which(ginfo$type %in% c("exon","5UTR","five_prime_utr","CDS","3UTR","three_prime_utr")),] %>%
              dplyr::group_by(gene_name,gene_id,transcript_id,type) %>%
              dplyr::summarise(typelen = sum(width)) %>%
              tidyr::spread(.,type,typelen,fill = 0) %>%
              dplyr::mutate(gtype = ifelse(`5UTR` > 0 | `CDS` > 0,"CD","NC")) %>%
              dplyr::mutate(cdsst = ifelse(`gtype` == "CD",`5UTR` + 1,1),
                            cdsed = ifelse(`gtype` == "CD",`cdsst` + `CDS`,exon),
                            tname = paste(gene_name,gene_id,transcript_id,cdsst,cdsed,exon,gtype,sep = "|"))

            # ===============================================================================
            # filter
            if(topN == 0){
              final.res <- leninfo
            }else{
              # choose type
              if(selecType == "lt"){
                final.res <- leninfo %>%
                  dplyr::group_by(gene_name,gene_id) %>%
                  dplyr::arrange(dplyr::desc(exon),dplyr::desc(CDS)) %>%
                  dplyr::slice_head(n = as.numeric(topN))
              }else if(selecType == "lcds"){
                final.res <- leninfo %>%
                  dplyr::group_by(gene_name,gene_id) %>%
                  dplyr::arrange(dplyr::desc(CDS),dplyr::desc(exon)) %>%
                  dplyr::slice_head(n = as.numeric(topN))
              }else{
                message("Please choose 'lt' or 'lcds'!")
              }
            }

            return(final.res)
          })



# setMethod("superExtract",
#           signature(object = "GenomeGTF"),
#           function(object,
#                    type = c("super.getLongName","super.filterRepTrans"),
#                    geneName = NULL,geneId = NULL,transId = NULL,
#                    selecType = c("lcds","lt"),topN = 1,sep = "|",
#                    ...){
#             type <- match.arg(type)
#
#             # select info
#             if(!is.null(object@representTrans)){
#               # load gtf
#               gtf <- object@representTrans
#
#               if(!is.null(geneName) & is.null(geneId) & is.null(transId)){
#                 ginfo <- gtf[which(gtf$gname %in% geneName),]
#               }else if(is.null(geneName) & !is.null(geneId) & is.null(transId)){
#                 ginfo <- gtf[which(gtf$gid %in% geneId),]
#               }else if(is.null(geneName) & is.null(geneId) & !is.null(transId)){
#                 ginfo <- gtf[which(gtf$tid %in% transId),]
#               }else if(!is.null(geneName) & is.null(geneId) & !is.null(transId)){
#                 ginfo <- gtf[which(gtf$gname == geneName & gtf$tid == transId),]
#               }else if(is.null(geneName) & !is.null(geneId) & !is.null(transId)){
#                 ginfo <- gtf[which(gtf$gid == geneId & gtf$tid == transId),]
#               }else{
#                 message("Please choose again.")
#               }
#
#               # check type
#               if(type == "super.getLongName"){
#                 if(sep != "|"){
#                   ginfo.sub <- ginfo[,c("gname","gid","tid","tname")] %>%
#                     dplyr::mutate(tname = gsub(pattern = "\\|",replacement = sep,x = .))
#                 }else{
#                   ginfo.sub <- ginfo[,c("gname","gid","tid","tname")]
#                 }
#                 return(ginfo.sub)
#               }else if(type == "super.filterRepTrans"){
#                 ginfo.sub <- ginfo[,c("gname","gid","tid","tlen","cdslen","tname")]
#
#                 # loop
#                 # plyr::ldply(unique(ginfo.sub$gid),function(x){
#                 #   tmp <- ginfo.sub[which(ginfo.sub$gid == x),]
#                 #
#                 #   # choose type
#                 #   if(selecType == "lt"){
#                 #     rankTran <- tmp %>%
#                 #       dplyr::arrange(dplyr::desc(tlen),dplyr::desc(cdslen))
#                 #   }else if(selecType == "lcds"){
#                 #     rankTran <- tmp %>%
#                 #       dplyr::arrange(dplyr::desc(cdslen),dplyr::desc(tlen))
#                 #   }else{
#                 #     message("Please choose 'lt' or 'lcds'!")
#                 #   }
#                 #
#                 #   # return length info
#                 #   if(topN == 0){
#                 #     return(rankTran)
#                 #   }else{
#                 #     rankTran <- rankTran %>%
#                 #       dplyr::slice_head(n = as.numeric(topN))
#                 #     return(rankTran)
#                 #   }
#                 # }) %>% data.frame() -> final.res
#
              #   if(topN == 0){
              #     final.res <- ginfo.sub
              #   }else{
              #     # choose type
              #     if(selecType == "lt"){
              #       final.res <- ginfo.sub %>%
              #         dplyr::group_by(gname,gid) %>%
              #         dplyr::arrange(dplyr::desc(tlen),dplyr::desc(cdslen)) %>%
              #         dplyr::slice_head(n = as.numeric(topN))
              #     }else if(selecType == "lcds"){
              #       final.res <- ginfo.sub %>%
              #         dplyr::group_by(gname,gid) %>%
              #         dplyr::arrange(dplyr::desc(cdslen),dplyr::desc(tlen)) %>%
              #         dplyr::slice_head(n = as.numeric(topN))
              #     }else{
              #       message("Please choose 'lt' or 'lcds'!")
              #     }
              #   }
              #
              #   return(final.res)
              # }else{
              #   message("Please choose correct type.")
              # }
#
#             }else{
#               message("representTrans slot is null.")
#             }
#           })

# getFeatureFromGenome
# extract "5UTR","five_prime_utr","3UTR","three_prime_utr","exon","CDS" and
# "intron" features sequnece from genome file.
#' Extract genomic features from GenomeGTF object
#'
#' This function extracts genomic features, such as 5' UTR, 3' UTR, exons, CDS and
#' introns, from a GenomeGTF object and returns the corresponding genomic
#' sequences. The genomic sequences are obtained by fetching the genome sequences
#' using the coordinates of the features.
#'
#' @param object A GenomeGTF object.
#' @param nameData A character vector or a data frame specifying the transcript
#' IDs to be extracted. If NULL, the function will use the geneName, geneId, and
#' transId arguments to filter the transcripts to be extracted. Default is NULL.
#' @param sep A character string specifying the separator used to concatenate
#' the gene name, gene ID, and transcript ID into a single name for the output
#' sequences. Default is "|".
#' @param geneName A character vector specifying the gene names to filter the
#' transcripts to be extracted. Default is NULL.
#' @param geneId A character vector specifying the gene IDs to filter the
#' transcripts to be extracted. Default is NULL.
#' @param transId A character vector specifying the transcript IDs to filter the
#' transcripts to be extracted. Default is NULL.
#' @param type A character string specifying the type of genomic features to be
#' extracted. Can be one of "5UTR", "five_prime_utr", "3UTR", "three_prime_utr",
#' "exon", "CDS", or "intron".
#' @param geneSeq A logical value specifying whether to extract the entire gene
#' sequence instead of the feature sequence. Default is FALSE.
#'
#' @importFrom plyr llply
#' @return A \code{DNAStringSet} object containing the genomic sequences of the
#' extracted features.
#' @method getFeatureFromGenome GenomeGTF
#' @export
setMethod("getFeatureFromGenome",
          signature(object = "GenomeGTF"),
          function(object,
                   nameData = NULL,
                   sep = "|",
                   geneName = NULL,geneId = NULL,transId = NULL,
                   type = c("5UTR","five_prime_utr","3UTR","three_prime_utr","exon","CDS","intron"),
                   geneSeq = FALSE){
            type <- match.arg(type)

            # load genome sequences
            genome <- object@genome

            # how to get id
            if(!is.null(nameData)){
              # load GTF
              if(type == "intron"){
                ginfo <- getIntronInfo(object = object,geneName = geneName,geneId = geneId,transId = nameData)
              }else{
                ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = nameData)
              }

              tid <- nameData
              tid.name <- names(tid)
            }else{
              # load GTF
              if(type == "intron"){
                ginfo <- getIntronInfo(object = object,geneName = geneName,geneId = geneId,transId = transId)
              }else{
                ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)
              }

              # get tids
              tid <- unique(ginfo$transcript_id)
            }

            # getseq function
            getMyseq <- function(genome,chr,start,end,strand){
              if(strand == "+"){
                seq <- genome[[chr]][start:end]
              }else{
                seq <- Biostrings::reverseComplement(genome[[chr]][start:end])
              }
              return(as.character(seq))
            }

            # filter chromosomes
            genomeChr <- names(mytest@genome)
            filteredInfo <- ginfo[which(ginfo$seqnames %in% genomeChr),]

            # get type info
            seqinfo.tmp <- filteredInfo[which(filteredInfo$type %in% c(type)),] %>%
              dplyr::select(gene_name,gene_id,transcript_id,type,seqnames,start,end,strand) %>%
              dplyr::group_by(gene_name,gene_id,transcript_id,type)

            # wthether get intron for other types
            if(geneSeq == TRUE){
              seqinfo <- seqinfo.tmp %>%
                dplyr::summarise(start = min(start),end = max(end),strand = unique(strand)) %>%
                dplyr::arrange(ifelse(strand == "+",c(start,end),c(desc(start),desc(end))),.by_group = T) %>%
                dplyr::mutate(typeNum = 1:dplyr::n())
            }else{
              seqinfo <- seqinfo.tmp %>%
                dplyr::arrange(ifelse(strand == "+",c(start,end),c(desc(start),desc(end))),.by_group = T) %>%
                dplyr::mutate(typeNum = 1:dplyr::n())
            }

            # tname
            all.res <- getTransInfo(object = object,geneName = unique(filteredInfo$gene_name),topN = 0) %>%
              dplyr::ungroup() %>%
              dplyr::select(transcript_id,tname)

            # merge
            seqinfo <- suppressMessages(dplyr::left_join(seqinfo,all.res))

            id <- unique(seqinfo$transcript_id)

            # progress bar
            pb <- progress::progress_bar$new(
              format = 'getFeatureFromGenome is running [:bar] :percent in :elapsed',
              total = length(id), clear = FALSE, width = 80
            )

            # loop get seq
            # x = 1
            purrr::map(1:length(id),function(x){
              pb$tick()
              tmp <- seqinfo[which(seqinfo$transcript_id == id[x]),]
              purrr::map(1:nrow(tmp),function(x){
                tmp1 <- tmp[x,]
                getMyseq(genome,tmp1$seqnames,tmp1$start,tmp1$end,tmp1$strand)
              }) -> seqlist

              if(type == "intron"){
                seq <- seqlist
              }else{
                seq <- paste0(seqlist,collapse = "")
              }

              Sys.sleep(0.05)
              return(seq)
            }) -> alllist

            # assign names
            if(type == "intron"){
              alllist <- unlist(alllist)
              names(alllist) <- paste(seqinfo$tname,seqinfo$typeNum,sep = sep)
            }else{
              names(alllist) <- unique(seqinfo$tname)
            }

            # to DNAStringSet
            res.seq <- Biostrings::DNAStringSet(unlist(alllist))
            return(res.seq)
          })


# getIntroInfo
# extract intron information from annotation file.
#' Get intron information from a GenomeGTF object
#'
#' This function extracts intron information from a GenomeGTF object for a
#' given gene or transcript.
#'
#' @param object A GenomeGTF object
#' @param geneName A character string specifying the gene name for which to
#' extract intron information. Default is NULL.
#' @param geneId A character string specifying the gene ID for which to extract
#'  intron information. Default is NULL.
#' @param transId A character string specifying the transcript ID for which to
#' extract intron information. Default is NULL.
#'
#' @return A data.frame containing the start and end coordinates of each intron,
#'  as well as other information such as the transcript ID and strand.
#' @method getIntronInfo GenomeGTF
#' @export
setMethod("getIntronInfo",
          signature(object = "GenomeGTF"),
          function(object,
                   geneName = NULL,geneId = NULL,transId = NULL){
            ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)

            # get tids
            tid <- unique(ginfo$transcript_id)

            # progress bar
            pb <- progress::progress_bar$new(
              format = 'getIntronInfo is running [:bar] :percent in :elapsed',
              total = length(tid), clear = FALSE, width = 80
            )

            # loop
            plyr::ldply(1:length(tid),function(x){
              pb$tick()

              tmp <- ginfo[which(ginfo$transcript_id == tid[x] & ginfo$type == "exon"),]

              # check whether only one exon
              if(nrow(tmp) == 1){
                message(paste0(tid[x]," has only one exon, skiping.",collapse = ""))
              }else{
                n.loop <- nrow(tmp)

                # get intron pos
                tmp.pos <- c(tmp$start[-1],tmp$end[-n.loop]) %>% sort()

                if(unique(tmp$strand) == "+"){
                  intron.info <- tmp[1:(n.loop - 1),]
                }else{
                  intron.info <- tmp[2:n.loop,]
                }

                # get index
                index <- c(1:length(tmp.pos))
                st.index <- index[index %% 2 == 1]
                sp.index <- st.index + 1

                # new data.frame for intron
                intron.info <- intron.info %>%
                  dplyr::mutate(start = tmp.pos[st.index] + 1,
                                end = tmp.pos[sp.index] - 1,
                                type = "intron") %>%
                  dplyr::mutate(width = end - start + 1)

                return(intron.info)
              }

              Sys.sleep(0.05)
            }) -> intronAll.info
            return(intronAll.info)
          })
