from Bio import SeqIO
from sys import argv

#Parses a -150 to +70 promoter sequences for each CDS
def user_seq(genbank_file):
    for index, record in enumerate(SeqIO.parse(open(genbank_file), "genbank")):
        DNA = record.seq
        for subfeatures in record.features:
            if subfeatures.type == 'CDS':
                try: gene_names = str(subfeatures.qualifiers['old_locus_tag'][0])
                except KeyError: gene_names = str(subfeatures.qualifiers['locus_tag'][0])
                orientation = subfeatures.strand
                upstream = int(subfeatures.location.nofuzzy_start) - 150
                downstream = int(subfeatures.location.nofuzzy_start) + 70
                sequences = DNA[upstream:downstream]
                print (">" + str(gene_names) + " " + str(orientation) + "\n" + str(sequences) + "\n"))
                    
user_seq(argv[1])
