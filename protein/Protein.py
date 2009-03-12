from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam, ProtParamData
from elm import ELM, Match
import re
import os
import time
import subprocess

class Protein:
	def __init__(self, record):
		self.sequence = record.seq
		self.id = record.id
		self.pa = ProtParam.ProteinAnalysis(str(self.sequence))
		self.features = {}
		self.scaled_features = {}
		self.findAminoAcidFeatures()
		self.findSecondaryStructureFeatures()
		self.findHydropathyFeatures()
		self.findInstabilityIndex()
		self.findELMFeatures()
		self.findPestFeatures()
		os.remove('temp.fa')
		os.remove('temp.ss')
		os.remove('temp.pest')
	
	def findELMFeatures(self):
		e = ELM.ELM()
		elms = e.elmList
		elmNames = elms.keys()
		for elmName in elmNames:
			elmMotif = elms[elmName]
			matcher = Match.Match(str(self.sequence), elmMotif)
			self.features[elmName] = matcher.getNumberMatches()


	def findInstabilityIndex(self):
		self.features['instability'] = self.pa.instability_index()

	def findAminoAcidFeatures(self):
		total = float(len(self.sequence))
		test_seq = str(self.sequence).upper()
		string = 'ACDEFGHIKLMNPQRSTUVWY'
		for char in string:
			self.features[char] = float(len(re.compile(char).findall(test_seq)))
		self.features['tiny'] = (self.features['A']+self.features['C']+self.features['U']+self.features['G']+self.features['S']+self.features['T'])/total
		self.features['small'] = (self.features['A']+self.features['C']+self.features['U']+self.features['D']+self.features['G']+self.features['N']+self.features['P']+self.features['S']+self.features['T']+self.features['V'])/total
		self.features['aliphatic'] = (self.features['A']+self.features['I']+self.features['L']+self.features['V'])/total
		self.features['aromatic'] = (self.features['F']+self.features['H']+self.features['W']+self.features['Y'])/total
		self.features['nonPolar'] = (self.features['A']+self.features['C']+self.features['U']+self.features['F']+self.features['G']+self.features['I']+self.features['L']+self.features['M']+self.features['P']+self.features['V']+self.features['W']+self.features['Y'])/total
		self.features['polar'] = (self.features['D']+self.features['E']+self.features['H']+self.features['K']+self.features['N']+self.features['Q']+self.features['R']+self.features['S']+self.features['T'])/total
		self.features['charged'] = (self.features['D']+self.features['E']+self.features['H']+self.features['K']+self.features['R'])/total
		self.features['basic'] = (self.features['H']+self.features['K']+self.features['R'])/total
		self.features['acidic'] = (self.features['D']+self.features['E'])/total
		self.features['length'] = total
		for char in string:
			self.features[char] = self.features[char]/total
		self.features['MW'] = self.pa.molecular_weight()

	def findHydropathyFeatures(self):
		self.features['pI'] = self.pa.isoelectric_point()
		hydropathy_total = 0.0
		hydropathy_plot = self.pa.protein_scale(ProtParamData.kd,10)
		for point in hydropathy_plot:
			hydropathy_total = hydropathy_total + point
		hydropathy_mean = hydropathy_total / float(len(hydropathy_plot))
		self.features['hydropathy'] = hydropathy_mean
	
	def findSecondaryStructureFeatures(self):
		fh = open('temp.fa', 'w')
		fh.write(">"+self.id+"\n")
		fh.write(str(self.sequence))
		fh.close()
		command = 'garnier -sequence temp.fa -outfile temp.ss > /dev/null 2>&1'
		sts = subprocess.Popen(command, shell=True)
		while os.path.exists('temp.ss') == False:
			time.sleep(0.1)
		oh = open('temp.ss', 'r')
		lines = oh.readlines()
		oh.close()
		p = re.compile('Residue totals: H:\s*(\d+)\s+E:\s*(\d+)')
		for line in lines:
			if p.search(line):
				m = p.search(line)
				helix = m.group(1)
				sheet = m.group(2)
		self.features['helix'] = float(helix) / float(len(str(self.sequence)))
		self.features['sheet'] = float(sheet) / float(len(str(self.sequence)))

	def findPestFeatures(self):
		command = "epestfind -sequence temp.fa -window 10 -outfile temp.pest -graph none -order score > /dev/null 2>&1"
		#sts = subprocess.call(command, shell=True)
		os.system(command)
		while os.path.exists('temp.pest') == False:
			time.sleep(0.1)
		oh = open('temp.pest', 'r')
		lines = oh.readlines()
		oh.close()
		important_line = lines[2]
		p1 = re.compile("^\s+(\d+)")
		p2 = re.compile("^\s+No")
		if p1.match(important_line):
			no_pest = p1.match(important_line).group(1)
		elif p2.match(important_line):
			no_pest = 0
		else:
			print 'eh?'
		self.features['pest'] = no_pest
#first digit here is number of PEST sites - if 0 will be 'No'

	
