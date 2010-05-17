from protein import Protein
from Bio import SeqIO
from svm import *
import random

def start():
    wt_proteins = getProteinData('data/wt.fa')
    eff_proteins = getProteinData('data/eff.fa')
    all_proteins = wt_proteins+eff_proteins
    keyList = all_proteins[0].features.keys()
    random.shuffle(all_proteins)
    #print all_proteins
    sortedKeys = sorted(keyList)
    excludedKeys = []
    for key in sortedKeys:
        max = float(all_proteins[0].features[key])
        min = float(all_proteins[0].features[key])
        for protein in all_proteins:
            if float(protein.features[key]) < min:
                min = float(protein.features[key])
            if float(protein.features[key]) > max:
                max = float(protein.features[key])
        for protein in all_proteins:
            x = float(protein.features[key])
            try:
                protein.scaled_features[key] = (x - min)/(max - min)
            except ZeroDivisionError:
                excludedKeys.append(key)
    labels = []
    samples = []
    for protein in all_proteins:
        label = 1
        for wt_protein in wt_proteins:
            if protein.id == wt_protein.id:
                label = 0
        labels.append(label)
        print label,
        sample = []
        i = 1
        for key in sortedKeys:
            exclude = 0
            for k in excludedKeys:
                if key == k:
                    exclude = 1
            if exclude == 0:
                sample.append(protein.scaled_features[key])
                print str(i)+":"+str(protein.scaled_features[key]),
                i += 1
        samples.append(sample)
        print
    
    problem = svm_problem(labels, samples);
    size = len(samples)
    kernels = [LINEAR, POLY, RBF]
    kname = ['linear','polynomial','rbf']
    param = svm_parameter(C = 10,nr_weight = 2,weight_label = [1,0],weight = [10,1])
    #for k in kernels:
    #    param.kernel_type = k;
    #    model = svm_model(problem,param)
    #    errors = 0
    #    for i in range(size):
    #        prediction = model.predict(samples[i])
    #        probability = model.predict_probability
    #        if (labels[i] != prediction):
    #            errors = errors + 1
    #    print "##########################################"
    #    print " kernel %s: error rate = %d / %d" % (kname[param.kernel_type], errors, size)
    #    print "##########################################"
    #param = svm_parameter(kernel_type = RBF, C=10)
    #model = svm_model(problem, param)
    #print "##########################################"
    #print " Decision values of predicting %s" % (samples[0])
    #print "##########################################"

    #print "Numer of Classes:", model.get_nr_class()
    #d = model.predict_values(samples[0])
    #for i in model.get_labels():
    #    for j in model.get_labels():
    #        if j>i:
    #            print "{%d, %d} = %9.5f" % (i, j, d[i,j])

    #param = svm_parameter(kernel_type = RBF, C=10, probability = 1)
    #model = svm_model(problem, param)
    #pred_label, pred_probability = model.predict_probability(samples[1])
    #print "##########################################"
    #print " Probability estimate of predicting %s" % (samples[1])
    #print "##########################################"
    #print "predicted class: %d" % (pred_label)
    #for i in model.get_labels():
    #    print "prob(label=%d) = %f" % (i, pred_probability[i])
    #model.save("test.model")

def getProteinData(file):
    fH = open(file, 'r')
    sequences = SeqIO.parse(fH, 'fasta')
    proteins = []
    for sequence in sequences:
        protein = Protein.Protein(sequence)
        proteins.append(protein)
    return proteins

if __name__ == "__main__":
    start()
