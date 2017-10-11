from sys import argv
from Bio import SeqIO

def get_16s(genbank_file, result_file):
  dna_file=open(result_file, "a")
  for index, record in enumerate(SeqIO.parse(open(genbank_file), "genbank")):
    DNA = record.seq
    for subfeatures in record.features:
      if subfeatures.type == "rRNA":
        product = subfeatures.qualifiers['product'] 
        split_product = str(product).split(" ")
        split_product = str(split_product).replace('"', "")
        split_product = str(split_product).replace("'", "")
        if "16S" in split_product:
          gene_name = genbank_file.replace(".gb", "_16S")
          start_coord = subfeatures.location.nofuzzy_start
          end_coord = subfeatures.location.nofuzzy_end
          feature_seq = DNA[start_coord:end_coord]
          if subfeatures.strand == -1:
            dna_seq = feature_seq.reverse_complement()
            dna_file.write(">" + str(gene_name) + "\n" + str(dna_seq) + '\n')
          else:
            dna_seq = feature_seq
            dna_file.write(">" + str(gene_name) + "\n" + str(dna_seq) + '\n')

get_16s(argv[1], "16s_rRNA.txt")