from Bio.Seq import Seq

class Protein:
	def __init__(self, record):
		self.sequence = record.seq
		self.id = record.id

