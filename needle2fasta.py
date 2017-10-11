from sys import argv
from Bio import AlignIO

fasta_file = argv[1].replace('.needle', '')

fasta_file = open(fasta_file+'.fasta', 'a')

for entries in AlignIO.parse(first_file, "emboss"):
	output = AlignIO.write(entries, results_files, 'fasta')

print "done"
		
