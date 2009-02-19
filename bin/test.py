from Bio import Application
from Bio.Application import _Option

class WaterCommandline(Application.AbstractCommandline):
	def __init__(self, cmd = "water"):
		Application.AbstractCommandline.__init__(self)
		self.program_name = cmd
		self.parameters = \
		[_Option(["-asequence"], ["input", "file"], None, 1,
		"First sequence to align"), 
		_Option(["-bsequence"], ["input", "file"], None, 1,
		"Second sequence to align"),
		_Option(["-gapopen"], ["input"], None, 1,
		"Gap open penalty"),
		_Option(["-gapextend"], ["input"], None, 1,
		"Gap extension penalty"),
		_Option(["-outfile"], ["output", "file"], None, 1,
		"Output file for the alignment"),
		_Option(["-datafile"], ["input", "file"], None, 0,
		"Matrix file"),
		_Option(["-similarity"], ["input"], None, 0,
		"Display percent identity and similarity"),
		_Option(["-nosimilarity"], ["input"], None, 0,
		"Do not display percent identity and similarity"),
		_Option(["-aformat"], ["input"], None, 0,
		"Display output in a different specified output format")] 

class PepstatsCommandline(Application.AbstractCommandline):
	def __init__(self, cmd = "pepstats"):
		Application.AbstractCommandline.__init__(self)
		self.program_name = cmd
		self.parameters = \
		[_Option(["-sequence"], ["input", "file"], None, 1,
		"First sequence to align"), 
		_Option(["-outfile"], ["output", "file"], None, 1,
		"Output file for the Statistics"),
		_Option(["-aadata"], ["input", "file"], None, 0,
		"Matrix file")]

w = WaterCommandline()
w.set_parameter("-asequence", "/home/sjcockell/P07830.fasta")
w.set_parameter("-bsequence", "/home/sjcockell/P07830.fasta")
w.set_parameter("-gapopen", "10.0")
w.set_parameter("-gapextend", "1.0")
w.set_parameter("-outfile", "out.txt")
Application.generic_run(w)
p = PepstatsCommandline()
p.set_parameter("-sequence", "/home/sjcockell/P07830.fasta")
p.set_parameter("-outfile", "stats.txt")
Application.generic_run(p)

