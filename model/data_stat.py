from Bio import SeqIO

def count_bases(fasta_file):
    bases_count = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    with open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequence = str(record.seq)
            for base in bases_count:
                bases_count[base] += sequence.count(base)

    return bases_count

def analyze_fasta(file_path):
    sequence_lengths = []

    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequence_length = len(record.seq)
            sequence_lengths.append(sequence_length)

    print(f"sequences number: {len(sequence_lengths)}")
    for i, length in enumerate(sequence_lengths, start=1):
        print(f"sequence {i} length is {length}。")

fasta_file_path = './data/human_2000.fasta'

result = count_bases(fasta_file_path)

print("A:", result['A'])
print("T:", result['T'])
print("C:", result['C'])
print("G:", result['G'])

analyze_fasta(fasta_file_path)
