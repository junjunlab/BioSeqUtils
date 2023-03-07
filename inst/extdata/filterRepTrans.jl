using GFF3
using DataStructures

function filterRepTrans(inputFile,outputFile,sep = "|")
    utr5Dict = Dict{String,Int64}()
    cdsDict = Dict{String,Int64}()
    exonDict = Dict{String,Int64}()

    # Open a GFF3 file.
    reader = open(GFF3.Reader, inputFile)

    ############################################################################################################
    # Iterate over GFF file
    ############################################################################################################
    for record in reader
        # save gene_id and gene name
        type = GFF3.featuretype(record)
        gene_name = split(GFF3.attributes(record,"gene_name")[1],";")[1]
        gene_id = GFF3.attributes(record,"gene_id")[1]
        transcript_id = GFF3.attributes(record,"transcript_id")[1]

        # chrom = GFF3.seqid(record)
        # strand = GFF3.strand(record)

        # dict_id = "$gene_name|$gene_id|$transcript_id"
        dict_id = join([gene_name,gene_id,transcript_id],sep)

        ######################################################################
        # calculate utr5 length
        if type in ["5UTR","five_prime_utr"]
            utr5Length = abs(GFF3.seqend(record) - GFF3.seqstart(record)) + 1
        else
            utr5Length = 0
        end
        # save in dict
        if !haskey(utr5Dict,dict_id)
            utr5Dict[dict_id] = utr5Length
        else
            utr5Dict[dict_id] += utr5Length
        end

        ######################################################################
        # calculate cds length
        if type in ["CDS"]
            CDSLength = abs(GFF3.seqend(record) - GFF3.seqstart(record)) + 1
        else
            CDSLength = 0
        end
        # save in dict
        if !haskey(cdsDict,dict_id)
            cdsDict[dict_id] = CDSLength
        else
            cdsDict[dict_id] += CDSLength
        end

        ######################################################################
        # calculate exon length
        if type in ["exon"]
            exonLength = abs(GFF3.seqend(record) - GFF3.seqstart(record)) + 1
        else
            exonLength = 0
        end
        # save in dict
        if !haskey(exonDict,dict_id)
            exonDict[dict_id] = exonLength
        else
            exonDict[dict_id] += exonLength
        end
    end


    ############################################################################################################
    # position dict
    ############################################################################################################
    posDict = Dict{String,Array{Int64}}()

    for (id,length) in utr5Dict
        if haskey(cdsDict,id)
            if haskey(exonDict,id)
                posDict[id] = [utr5Dict[id] + 1,utr5Dict[id] + cdsDict[id],exonDict[id]]
            end
        end
    end


    outputfile = open(outputFile,"w")

    for (id,length) in posDict
        gname,gid,tid = split(id,sep)
        cdsst,cdsed,exonlen = length
        cdslen = cdsDict[id]
        write(outputfile,join([gname,gid,tid,cdsst,cdsed,exonlen,id,cdslen],"\t")*"\n")
    end

    # close file
    close(outputfile)
end