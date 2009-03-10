from Bio.Seq import Seq
import re
import os
import time
import subprocess

class Protein:
	def __init__(self, record):
		self.sequence = record.seq
		self.id = record.id
		self.features = {}
	
	def findELMFeatures(self):
		pass
	
	def findAminoAcidFeatures(self):
		total = float(len(self.sequence))
		test_seq = str(self.sequence).upper()
		ala = float(len(re.compile('A').findall(test_seq)))
		asx = float(len(re.compile('B').findall(test_seq)))
		cys = float(len(re.compile('C').findall(test_seq)))
		asp = float(len(re.compile('D').findall(test_seq)))
		glu = float(len(re.compile('E').findall(test_seq)))
		phe = float(len(re.compile('F').findall(test_seq)))
		gly = float(len(re.compile('G').findall(test_seq)))
		his = float(len(re.compile('H').findall(test_seq)))
		ile = float(len(re.compile('I').findall(test_seq)))
		lys = float(len(re.compile('K').findall(test_seq)))
		leu = float(len(re.compile('L').findall(test_seq)))
		met = float(len(re.compile('M').findall(test_seq)))
		asn = float(len(re.compile('N').findall(test_seq)))
		pro = float(len(re.compile('P').findall(test_seq)))
		gln = float(len(re.compile('Q').findall(test_seq)))
		arg = float(len(re.compile('R').findall(test_seq)))
		ser = float(len(re.compile('S').findall(test_seq)))
		thr = float(len(re.compile('T').findall(test_seq)))
		sec = float(len(re.compile('U').findall(test_seq)))
		val = float(len(re.compile('V').findall(test_seq)))
		trp = float(len(re.compile('W').findall(test_seq)))
		xaa = float(len(re.compile('X').findall(test_seq)))
		tyr = float(len(re.compile('Y').findall(test_seq)))
		glx = float(len(re.compile('Z').findall(test_seq)))
		self.features['A'] = ala/total
		self.features['C'] = cys/total
		self.features['D'] = asp/total
		self.features['E'] = glu/total
		self.features['F'] = phe/total
		self.features['G'] = gly/total
		self.features['H'] = his/total
		self.features['I'] = ile/total
		self.features['K'] = lys/total
		self.features['L'] = leu/total
		self.features['M'] = met/total
		self.features['N'] = asn/total
		self.features['P'] = pro/total
		self.features['Q'] = gln/total
		self.features['R'] = arg/total
		self.features['S'] = ser/total
		self.features['T'] = thr/total
		self.features['U'] = sec/total
		self.features['V'] = val/total
		self.features['W'] = trp/total
		self.features['Y'] = tyr/total
		self.features['tiny'] = (ala+cys+sec+gly+ser+thr)/total
		self.features['small'] = (ala+cys+sec+asp+gly+asn+pro+ser+thr+val)/total
		self.features['aliphatic'] = (ala+ile+leu+val)/total
		self.features['aromatic'] = (phe+his+trp+tyr)/total
		self.features['nonPolar'] = (ala+cys+sec+phe+gly+ile+leu+met+pro+val+trp+tyr)/total
		self.features['polar'] = (asp+glu+his+lys+asn+gln+arg+ser+thr)/total
		self.features['charged'] = (asp+glu+his+lys+arg)/total
		self.features['basic'] = (his+lys+arg)/total
		self.features['acidic'] = (asp+glu)/total
		self.features['length'] = total

	def findHydropathyFeatures(self):
		pass
	
	def findSecondaryStructureFeatures(self):
		fh = open('temp.fa', 'w')
		fh.write(">"+self.id+"\n")
		fh.write(str(self.sequence))
		fh.close()
		command = 'garnier -sequence temp.fa -outfile temp.txt &> /dev/null'
		sts = subprocess.call(command, shell=True)
		while os.path.exists('temp.txt') == False:
			time.sleep(0.1)
		oh = open('temp.txt', 'r')
		lines = oh.readlines()
		oh.close()
		p = re.compile('Residue totals: H: (\d+)\s+E: (\d+)')
		for line in lines:
			if p.search(line):
				m = p.search(line)
				helix = m.group(1)
				sheet = m.group(2)
		self.features['helix'] = float(helix) / float(len(str(self.sequence)))
		self.features['sheet'] = float(sheet) / float(len(str(self.sequence)))
		os.remove('temp.fa')
		os.remove('temp.txt')

	def findPestFeatures(self):
		pass
	
