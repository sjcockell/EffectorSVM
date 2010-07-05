from protein import Protein
from Bio import SeqIO
from svm import *
import random

def main():
    model = svm_model('train_features_scaled.model')
    proteins = getProteinData('data/effectors_yersinia_pestis.fa')
    keyList = proteins[0].features.keys()
    random.shuffle(proteins)
    sortedKeys = sorted(keyList)
    excludedKeys = []
    samples = []
    for key in sortedKeys:
        max = float(proteins[0].features[key])
        min = float(proteins[0].features[key])
        for protein in proteins:
            if float(protein.features[key]) < min:
                min = float(protein.features[key])
            if float(protein.features[key]) > max:
                max = float(protein.features[key])
        for protein in proteins:
            x = float(protein.features[key])
            try:
                protein.scaled_features[key] = (x - min)/(max - min)
            except ZeroDivisionError:
                excludedKeys.append(key)
    for protein in proteins:
        sample = []
        i = 1
        for key in sortedKeys:
            exclude = 0
            for k in excludedKeys:
                if key == k:
                    exclude = 1
            if exclude == 0:
                sample.append(protein.scaled_features[key])
                i += 1
        samples.append(sample)
    #pred_label, pred_probability = model.predict_probability(proteins[1])
    #print "##########################################"
    #print " Probability estimate of predicting %s" % (proteins[1])
    #print "##########################################"
    #print "predicted class: %d" % (pred_label)
    #for i in model.get_labels():
    #    print "prob(label=%d) = %f" % (i, pred_probability[i])
    #for sample in samples:
    #    print model.predict(sample)
    print "Numer of Classes:", model.get_nr_class()
    d = model.predict_values(samples[0])
    for i in model.get_labels():
        for j in model.get_labels():
            if j>i:
                print "{%d, %d} = %9.5f" % (i, j, d[i,j])


def getProteinData(file):
    fH = open(file, 'r')
    sequences = SeqIO.parse(fH, 'fasta')
    proteins = []
    for sequence in sequences:
        protein = Protein.Protein(sequence)
        proteins.append(protein)
    return proteins

if __name__ == '__main__':
    main()
