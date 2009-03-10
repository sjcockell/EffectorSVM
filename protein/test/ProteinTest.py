from protein import Protein as pro
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os

record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
                    IUPAC.protein),
                    id="YP_025292.1", name="HokC",
                    description="toxic membrane protein, small")

protein = pro.Protein(record)
protein.findAminoAcidFeatures()
protein.findSecondaryStructureFeatures()
feature = protein.features
print feature
