from Bio import SeqIO
import random

fasta_file = "./data/xx.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))

random.shuffle(records)

output_file = "./data/xx_shuffle.fasta"
with open(output_file, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

print("Sequences shuffled and saved to:", output_file)