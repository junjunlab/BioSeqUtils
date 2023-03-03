globalVariables(c("cdslen", "end", "start", 'tga', "tlen"))

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


#' filterRepTrans method for GenomeGTF objects
#'
#' @title filterRepTrans method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("filterRepTrans",function(object,...) standardGeneric("filterRepTrans"))


#' getLongName method for GenomeGTF objects
#'
#' @title getLongName method for GenomeGTF objects
#' @param object A GenomeGTF object
#' @param ... Other arguments.
#' @export
setGeneric("getLongName",function(object,...) standardGeneric("getLongName"))


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

# ==============================================================================
# setMethod
# ==============================================================================
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

            cat("## representTrans slot is NULL.\n")
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


# filterRepTrans
# filter the longest transcript from annotation file and return transcript length
# and CDS length with data frame.
#' Filter redundant transcripts and retain the longest ones based on gene ID or
#' gene name
#'
#' This function filters redundant transcripts based on gene ID or gene name,
#' and retains the longest ones according to the selected criteria.
#'
#' @param object A GenomeGTF object or a data.frame with GTF
#'  format.
#' @param geneName A character vector of gene names. Default is NULL.
#' @param geneId A character vector of gene IDs. Default is NULL.
#' @param selecType A character vector specifying the selection criteria.
#' Must be either "lt" for longest transcript or "lcds" for longest CDS. Default
#' is "lcds".
#' @param topN An integer indicating the number of longest transcripts to retain.
#' Default is 1.
#' Set to 0 to retain all transcripts.
#' @param sep A character specifying the separator used to separate gene name and
#' transcript name.
#' Default is "|".
#' @return A data.frame containing the longest transcripts for each gene based on
#' the selected criteria.
#' @method filterRepTrans GenomeGTF
#' @export
setMethod("filterRepTrans",
          signature(object = "GenomeGTF"),
          function(object,
                   geneName = NULL,
                   geneId = NULL,
                   selecType = c("lcds","lt"),
                   # 0 for all selections
                   topN = 1,sep = "|"){
            # doing now
            selecType <- match.arg(selecType)

            # get data
            if(as.character(class(object)) == "GenomeGTF"){
              gtf <- as.data.frame(object@gtf)

              # choose type
              if(!is.null(geneName) & is.null(geneId)){
                tga <- gtf[which(gtf$gene_name %in% geneName),]
              }else if(is.null(geneName) & !is.null(geneId)){
                tga <- gtf[which(gtf$gene_id %in% geneId),]
              }else if(!is.null(geneName) & !is.null(geneId)){
                tga <- gtf[which(gtf$gene_name == geneName & gtf$gene_id == geneId),]
              }

            }else if(as.character(class(object)) == "data.frame"){
              tga <- object[which(object$gene_id %in% geneId),]
            }

            gid = unique(tga$gene_id)

            # progress bar
            pb <- progress::progress_bar$new(
              format = 'filterRepTrans is running [:bar] :percent in :elapsed',
              total = length(gid), clear = FALSE, width = 80
            )

            # loop
            plyr::ldply(1:length(gid),function(x){
              pb$tick()

              tg <- tga[which(tga$gene_id == gid[x]),]

              # calculate transcript/CDS length
              tida <- unique(tg$transcript_id)
              plyr::ldply(1:length(tida),function(x){

                tmp.info <- tg[which(tg$transcript_id == tida[x]),]
                t.type <- unique(tmp.info$type)

                # exon info
                tmp <- tmp.info[which(tmp.info$type == "exon"),]

                # check transcript whether has CDS
                if("CDS" %in% t.type){
                  tmp2 <- tmp.info[which(tmp.info$type == "CDS"),]
                  cdslen = sum(tmp2$width)
                }else{
                  cdslen = 0
                }

                # return data.frame
                tranLength <- data.frame(gid = unique(tg$gene_id),
                                         tid = tida[x],
                                         tlen = sum(tmp$width),
                                         cdslen = cdslen)

                return(tranLength)
              }) -> lenInfo

              # choose type
              if(selecType == "lt"){
                rankTran <- lenInfo %>%
                  dplyr::arrange(dplyr::desc(tlen),dplyr::desc(cdslen))
              }else if(selecType == "lcds"){
                rankTran <- lenInfo %>%
                  dplyr::arrange(dplyr::desc(cdslen),dplyr::desc(tlen))
              }else{
                message("Please choose 'lt' or 'lcds'!")
              }

              # get long name
              rankTran.name <- getLongName(object = object,transId = rankTran$tid,sep = sep)
              rankTran$tname <- rankTran.name$tname

              # return length info
              if(topN == 0){
                return(rankTran)
              }else{
                plyr::ldply(unique(rankTran$gid),function(x){
                  rankTran <- rankTran[which(rankTran$gid == x),] %>%
                    dplyr::slice_head(n = as.numeric(topN))
                  return(rankTran)
                }) -> rankTran
                return(rankTran)
              }

              Sys.sleep(0.05)
            }) -> rep.info

            return(rep.info)
          })


# getLongName
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
#' @param sep Separator character to use in the long name. Default is "|".
#'
#' @return A data frame with columns gid, tid, and tname.
#'
#' @importFrom plyr llply ldply
#' @importFrom dplyr mutate desc select arrange slice_head
#' @method getLongName GenomeGTF
#' @export
setMethod("getLongName",
          signature(object = "GenomeGTF"),
          function(object,
                   geneName = NULL,geneId = NULL,transId = NULL,
                   sep = "|"){
            # load GTF
            ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)

            # get tids
            tid <- unique(ginfo$transcript_id)

            # loop
            plyr::ldply(1:length(tid),function(x){
              ginfo.s <- ginfo[which(ginfo$transcript_id == tid[x]),]
              gname <- ginfo.s$gene_name[1]
              tid <- ginfo.s$transcript_id[1]

              # transcript length
              tlen <- sum(ginfo.s[which(ginfo.s$type == "exon"),"width"])

              # check whether exist 5UTR and get CDS start position
              if("5UTR" %in% ginfo.s$type){
                cds.st <- sum(ginfo.s[which(ginfo.s$type == "5UTR"),"width"]) + 1
              }else{
                cds.st <- 1
              }

              # get CDS length
              if("CDS" %in% ginfo$type){
                cds.len <- sum(ginfo.s[which(ginfo.s$type == "CDS"),"width"]) + 1
              }else{
                cds.len <- tlen
              }

              # get CDS stop position
              cds.sp <- cds.len + cds.st - 1

              # name
              tname <- paste(gname,tid,cds.st,cds.sp,tlen,sep = sep)
              res <- data.frame(gid = gname,tid = tid,tname = tname)
              return(res)
            }) -> resName
            return(resName)
          })


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
            myFASTA <- object@genome

            # how to get id
            if(!is.null(nameData)){
              # load GTF
              ginfo <- as.data.frame(object@gtf)

              # character or data.frame
              if(is.character(nameData)){
                tid <- nameData
                tid.name <- names(tid)
              }else{
                tid <- nameData$tid
                tid.name <- nameData$tname
              }
            }else{
              # load GTF
              ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)

              # get tids
              tid <- unique(ginfo$transcript_id)
            }

            # progress bar
            pb <- progress::progress_bar$new(
              format = 'getFeatureFromGenome is running [:bar] :percent in :elapsed',
              total = length(tid), clear = FALSE, width = 80
            )

            # loop
            plyr::llply(1:length(tid),function(x){
              pb$tick()

              # check utr or CDS
              if(type == "intron"){
                ginfo.fet <- getIntronInfo(object,transId = tid[x])
              }else{
                ginfo.fet <- ginfo[which(ginfo$transcript_id == tid[x]),]
              }

              # check features
              if(type %in% unique(ginfo.fet$type)){

                if(type == "intron"){
                  exon.info <- ginfo.fet[which(ginfo.fet$type %in% type),]
                }else{
                  exon.info <- ginfo[which(ginfo$type %in% type & ginfo$transcript_id == tid[x]),]
                }

                # check whether genome has this chromosome
                if(unique(exon.info$seqnames) %in% names(myFASTA)){
                  # strand
                  strand.info <- unique(exon.info$strand)
                  if(strand.info == "+"){
                    exon.arrange <- exon.info %>% dplyr::arrange(start,end)
                  }else{
                    exon.arrange <- exon.info %>% dplyr::arrange(dplyr::desc(start),dplyr::desc(end))
                  }

                  # whether extract gene sequnce
                  if(geneSeq == TRUE){
                    exon.arrange <- exon.arrange[1,] %>%
                      dplyr::mutate(start = exon.arrange$start[1],
                                    end = exon.arrange$end[nrow(exon.arrange)])
                  }else{
                    exon.arrange <- exon.arrange
                  }

                  # fetch exon seq
                  # x = 1
                  plyr::llply(1:nrow(exon.arrange),function(x){
                    tmp <- exon.arrange[x,]

                    if(strand.info == "+"){
                      seq <- myFASTA[[tmp$seqnames]][tmp$start:tmp$end]
                    }else{
                      seq <- Biostrings::reverseComplement(myFASTA[[tmp$seqnames]][tmp$start:tmp$end])
                    }

                    # as.character(seq)
                    as(seq,"character")
                  }) -> exon.seqs

                  # final.seq <- Biostrings::DNAStringSet(paste0(exon.seqs,collapse = ""))
                  if(is.null(nameData)){
                    if(type == "intron"){
                      final.seq <- exon.seqs

                      tid.name <- paste(unique(exon.arrange$gene_name),unique(exon.arrange$gene_id),
                                        unique(exon.arrange$transcript_id),1:nrow(exon.arrange),
                                        sep = sep)
                    }else{
                      final.seq <- paste0(exon.seqs,collapse = "")

                      tid.name <- paste(unique(exon.arrange$gene_name),unique(exon.arrange$gene_id),
                                        unique(exon.arrange$transcript_id),
                                        sep = sep)
                    }
                    # assign name
                    names(final.seq) <- tid.name
                  }else{
                    if(type == "intron"){
                      final.seq <- exon.seqs

                      # assign name
                      names(final.seq) <- paste(tid.name[x],exon.arrange$exon_number,sep = sep)
                    }else{
                      final.seq <- paste0(exon.seqs,collapse = "")

                      # assign name
                      names(final.seq) <- tid.name[x]
                    }
                  }

                  return(final.seq)
                }else{
                  message(paste0("Gene: ",tid[x]," chromosome ",unique(exon.info$seqnames),
                                 " is not in genome file. Skiping this."))
                }
              }else{
                message(paste0("Gene: ",tid[x]," ",type," is not found. Skiping this."))
              }

              # return(final.seq)
            }) -> all.seqs

            res.seq <- Biostrings::DNAStringSet(unlist(all.seqs))
            Sys.sleep(0.05)
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
