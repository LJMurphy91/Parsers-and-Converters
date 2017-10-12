from sys import argv
from Bio import SeqIO
import re

#USE - python annotated_inteins.py <genbank.gb>

def get_inteins(genbank_file, result_file):
  filename = genbank_file
  extension = re.compile("\.(.*)")
  strain = re.sub(extension, "", filename)
  outfile=open(result_file, "a")
  for index, record in enumerate(SeqIO.parse(open(genbank_file), "genbank")):
    DNA = record.seq
    for subfeatures in record.features:
      if subfeatures.type == "CDS":
        try: notes = subfeatures.qualifiers['note'][0]
        except KeyError: notes="no note"
        try: product = subfeatures.qualifiers["product"][0]
        except KeyError: product="no product"
        str_notes = "\t".join(notes)
        str_product = "\t".join(product)
        if "intein" in str_product:
          locus = subfeatures.qualifiers['db_xref'][0]
          try: proteins = subfeatures.qualifiers['translation'][0]
          except KeyError: proteins = DNA[subfeatures.location.nofuzzy_start:subfeatures.location.nofuzzy_end].translate()
          protein_seq = proteins
          outfile.write(">" + str(strain) + "_" + str(locus) + "\n" + str(protein_seq) + '\n')
        if "intein" in str_notes:
          locus = subfeatures.qualifiers['locus_tag'][0]
          try: proteins = subfeatures.qualifiers['translation'][0]
          except KeyError: proteins = DNA[subfeatures.location.nofuzzy_start:subfeatures.location.nofuzzy_end].translate()
          outfile.write(">" + str(strain) + "_" + str(locus) + "\n" + str(proteins) + '\n')

get_inteins(argv[1], "annotated_inteins.txt")