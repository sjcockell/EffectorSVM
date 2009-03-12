from Bio import SeqIO
import random
def main():
	fh1 = open('/home/sjcockell/git/elmBac/data/ecoli.fa', 'r')
	fh2 = open('/home/sjcockell/git/elmBac/data/epec.fa', 'r')
	fh3 = open('/home/sjcockell/git/elmBac/data/salmonella.fa', 'r')
	fh4 = open('/home/sjcockell/git/elmBac/data/shigella.fa', 'r')
	
	sequences1 = SeqIO.parse(fh1, 'fasta')
	sequences2 = SeqIO.parse(fh2, 'fasta')
	sequences3 = SeqIO.parse(fh3, 'fasta')
	sequences4 = SeqIO.parse(fh4, 'fasta')
	
	wt_sequences = {}
	effector_sequences = {}
	
	for sequence in sequences1:
		wt_sequences[sequence.id] = str(sequence.seq)
	for sequence in sequences2:
		effector_sequences[sequence.id] = str(sequence.seq)
	for sequence in sequences3:
		effector_sequences[sequence.id] = str(sequence.seq)
	for sequence in sequences4:
		effector_sequences[sequence.id] = str(sequence.seq)

	fh5 = open('data/wt.fa', 'w')
	wt_keys = wt_sequences.keys()
	keys = random.sample(wt_keys, 50)
	for key in keys:
		fh5.write(">"+key+"\n")
		fh5.write(wt_sequences[key]+"\n")

	fh6 = open('data/eff.fa', 'w')
	effector_keys = effector_sequences.keys()
	keys = random.sample(effector_keys, 50)
	for key in keys:
		fh6.write(">"+key+"\n")
		fh6.write(effector_sequences[key]+"\n")
	
	fh1.close()
	fh2.close()
	fh3.close()
	fh4.close()
	fh5.close()
	fh6.close()

if __name__ == "__main__":
	main()
