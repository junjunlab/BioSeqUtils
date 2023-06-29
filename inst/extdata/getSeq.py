from pyfaidx import Fasta
import re

# save id in dictionary
def id_dict(transcrip_id,new_id):
    idDict = {}
    i = 0
    for tid in transcrip_id:
        idDict[tid] = ">"+new_id[i]
        i += 1
    return idDict

# extract record from attributes
def extract_record(attribute_string):
    # gene_name_pattern = r'gene_name\s+"([^"]+)"'
    # gene_id_pattern = r'gene_id\s+"([^"]+)"'
    transcript_id_pattern = r'transcript_id\s+"([^"]+)"'
    
    # gene_name = re.search(gene_name_pattern, attribute_string)
    # gene_id = re.search(gene_id_pattern, attribute_string)
    transcript_id = re.search(transcript_id_pattern, attribute_string)
    
    # gene_name = gene_name.group(1) if gene_name else None
    # gene_id = gene_id.group(1) if gene_id else None
    transcript_id = transcript_id.group(1) if transcript_id else None
    
    return transcript_id

# extract sequence from genome
def extract_sequnence(gtfFile,genomeObj,target_idDict,type = "exon"):
    seq_dict = {value: "" for key, value in target_idDict.items()}

    with open(gtfFile) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            # attributes for gtf
            attributes = fields[8]
            # get feature type
            feature = fields[2]
            # get type of line
            if feature != type:
                continue
            # transcript id
            transcript_id = extract_record(attributes)
            # seq_dict key
            newkey = target_idDict.get(transcript_id)
            # get sequncence form genome sequence
            if newkey != None:
                chrom,start,end = fields[0],int(fields[3]),int(fields[4])
                
                # check whether genome has this chromosome
                if chrom in genomeObj.keys():
                    type_seq = genomeObj[chrom][start-1:end]
                    # wthether + strand or - strand
                    if fields[6] == "+":
                        seq_dict[newkey] += str(type_seq)
                    else:
                        type_seq = type_seq.complement.reverse
                        seq_dict[newkey] = str(type_seq) + str(seq_dict[newkey])
                else:
                    continue
    return seq_dict
                    

# output fasta file
def write_fasta(faDict,fileName,lineBases = 80):
    # output save in another file
    output_fa = open(fileName,'w')

    # separate sequences
    for key,val in faDict.items():
        if len(val) != 0:
            output_fa.write(key + '\n')
            while len(val) > lineBases:
                output_fa.write(val[0:lineBases] + '\n')
                val = val[lineBases:len(val)]
            output_fa.write(val + '\n')

    # file close
    output_fa.close()
                    
# one step to perform sequnence extraction from genome data
def py_extractSequence(gtf_file,genome_file,transcript_id,new_id,type,out_file):
    target_idDict = id_dict(transcript_id,new_id)
    genome_my = Fasta(genome_file, sequence_always_upper=True)
    seq_dict = extract_sequnence(gtf_file,genome_my,target_idDict,type=type)
    write_fasta(seq_dict,out_file)