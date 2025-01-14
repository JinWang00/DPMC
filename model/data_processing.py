import random
from Bio import SeqIO
from Bio.Seq import Seq

filename = './data/human_2000.fasta'
reads = []

with open(filename, 'r') as handle:
	for record in SeqIO.parse(handle, 'fasta'):
		reads.append(record)
print(len(reads))

data = list(range(len(reads)))
random.shuffle(data)
train = data[0:1500]
valid = data[1500:2000]

train_record = []
valid_record = []

def check(read):
	record = read
	my_dna = str(read.seq.upper()) # Convert lowercase to uppercase
	for i, base in enumerate(my_dna):
		if base not in 'ACGT':
			my_dna = my_dna.replace(base,'A') # Replaced invalid base with 'A'

	record.seq = Seq(my_dna)
	for i, base in enumerate(record.seq):
		if base not in 'ACGT':
			print(record.seq[i])
	return record
for i in train:
	read = check(reads[i])
	train_record.append(read)
for i in valid:
	read = check(reads[i])
	valid_record.append(read)

#save
SeqIO.write(train_record, "./data/human_1-1500.fasta", "fasta")
SeqIO.write(valid_record, "./data/human_1500-2000.fasta", "fasta")
print('finish')
