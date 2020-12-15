#setting up actual dataset sequences for Fragoza to embed
import pandas as pd
import sys
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
import os.path
from os import path
import itertools
import pickle
import copy
import re
from Bio.SeqIO import *
from Bio import SeqIO
import matplotlib.pyplot as plt

def openFasta(fastaName):
    record = SeqIO.read(fastaName, "fasta")
    return str(record.seq)


def mutateSequence(fastaName, mutation):
    # opens FASTA and mutates sequence, saves mutated Fasta in same folder
    # open FASTA
    record = SeqIO.read(fastaName, "fasta")
    sequenceToMutate = str(record.seq)
    halves = mutation.split(",")
    actualMutation = halves[1]
    firstAA = actualMutation[0]
    mutAA = actualMutation[len(actualMutation) - 1]
    number = int(actualMutation[1:len(actualMutation) - 1])
    # check that the AA at that position is correct
    # print ("len seq: ", len(sequenceToMutate))
    # print ("actualMutation: ", actualMutation)
    # print ("number: ", number)
    if number - 1 <= len(sequenceToMutate):
        # print ("at pos: ", sequenceToMutate[number - 1])
        if firstAA == sequenceToMutate[number - 1]:
            # print ("mutating")
            index = number - 1
            newSeq = sequenceToMutate[:index] + mutAA + sequenceToMutate[index + 1:]
            return newSeq
        else:
            return -1
    else:
        return -1

def flipDisruption(disr):
    if disr == 1:
        return 0
    else:
        return 1

def getFragozaInteractionDFs():
    #returns csv of interactions and sequences below <= 2000 AA
    # attaching NP to each row
    origDatasetDtypes = {
        "Chrom": str,
        "Pos": float,
        "dbSNP_id": str,
        "Mutation": str,
        "AF_all": float,
        "Target Entrez GeneID": str,
        "HGVS_cDNA": str,
        "UniProt": str,
        "Interactor Entrez GeneID": str,
        "Disruption": float,
    }
    # print (origDatasetDtypes)
    origDataset = pd.read_csv("41467_2019_11959_MOESM6_ESM.csv", sep=',')
    origDataset = origDataset.astype(origDatasetDtypes)

    origDataset['transcript'] = origDataset.apply(lambda row: row['HGVS_cDNA'].split(':')[0], axis=1)
    origDataset['mutationAA'] = origDataset.apply(lambda row: row['UniProt'].split(',')[1], axis=1)

    openedKey = pd.read_csv("transcriptProteinTableFragoza.csv")
    merged = pd.merge(origDataset, openedKey, how='left', on='transcript')
    merged['saveIDTarget'] = merged['Target Entrez GeneID'] + "_" + merged['mutationAA']
    print ("intermed length: ", merged.shape) #(4109, 15)
    seqTable = {"ID": [], "Seq": [], 'Len': []}
    # saveIDTarget (mutation partner) mut id
    # "Interactor Entrez GeneID" other
    print(merged.iloc[0])
    badMutants = []
    for i, row in merged.iterrows():
        # open orig FASTA
        origID = row['Target Entrez GeneID']
        mutation = row['mutationAA']
        np = row['product']
        fastaFolder = './FragozaFastas/'
        mutationU = row['UniProt']
        saveID = row['saveIDTarget']
        fastaName = fastaFolder + np + ".fa"
        newSeq = mutateSequence(fastaName, mutationU)
        wtSeq = openFasta(fastaName)
        if not newSeq == -1:
            seqTable["ID"].append(saveID)
            seqTable["Seq"].append(newSeq)
            seqTable['Len'].append(len(newSeq))
        else:
            #print("Bad mutation")
            badMutants.append(saveID)
        if len(wtSeq) > 1:
            seqTable["ID"].append(origID)
            seqTable["Seq"].append(wtSeq)
            seqTable['Len'].append(len(wtSeq))
    badNamesInteractors = [122183, 441521, 1409, 51207, 9465]

    # saving those interactor wt to a dataframe
    uniprotInteractors = pd.read_csv("withIDsInteractorUniprot.tab", sep='\t')
    # now for each not in the bad set, get the seqT
    for i, row in uniprotInteractors.iterrows():
        idInt = row['yourlist']
        sequence = row['Sequence']
        if idInt not in badNamesInteractors:
            seqTable['ID'].append(idInt)
            seqTable['Seq'].append(sequence)
            seqTable['Len'].append(len(sequence))

    badInteractors = ['390535', '541471', '729862', '202459', '155060', '401508', '80028']
    goodInteractors = [10597, 285733, 145946, 440184, 440321, 9753, 122183, 441521, 1409, 51207, 9465]

    for interactors in goodInteractors:
        fastaName = "./FragozaInteractors/" + str(interactors) + ".fasta"
        fastaSeq = openFasta(fastaName)
        seqTable['Seq'].append(fastaSeq)
        seqTable['ID'].append(interactors)
        seqTable['Len'].append(len(fastaSeq))

    asDF = pd.DataFrame(seqTable)
    print(asDF.drop_duplicates())

    # save as CSV
    # trim lengths here
    #print(asDF.shape)
    #print (asDF)
    #distribution of sequence lengths:
    asDF = asDF.drop_duplicates()
    print ("total sequences: ", asDF.shape)
    asDF = asDF[asDF['Len'] <= 2000]
    asDF = asDF.drop_duplicates()
    print ("below 2000 AA: ", asDF.shape)

    #total sequences:  (8927, 3)
    #below 2000 AA:  (8898, 3)

    print(asDF.shape)
    passedIDS = asDF['ID'].to_list()
    hist = asDF.hist()
    plt.show()
    asDF = asDF.drop(['Len'], axis=1)



    asDF.to_csv("FragozaSequencesMaxLen2000.tab", sep="\t", header=False, index=False)

    PIPRPairs = {'v1': [], 'v2': [], 'seq1':[], 'seq2':[], 'label': [], 'type': []}  # 1 if interacts, 0 else

    for i, row in merged.iterrows():
        # open orig FASTA
        print (row)
        origID = row['Target Entrez GeneID']  # wt 1
        saveID = row['saveIDTarget']  # mt 1
        wtID2 = row["Interactor Entrez GeneID"]
        if wtID2 in passedIDS and origID in passedIDS and saveID in passedIDS:
            disruption = row['Disruption']  # 0 if nondisrupting, 1 if disrupting
            PIPRPairs['v1'].append(origID)
            #get v1 seq
            #seq1 = asDF.loc[asDF.ID == origID]
            PIPRPairs['v2'].append(wtID2)
            PIPRPairs['label'].append(1)
            PIPRPairs['type'].append("WT")
            PIPRPairs['v1'].append(saveID)
            PIPRPairs['v2'].append(wtID2)
            #print("DIS WAS: ", disruption)
            #print("DIS NOW: ", flipDisruption(disruption))
            PIPRPairs['label'].append(flipDisruption(disruption))
            PIPRPairs['type'].append("MT")


    # 868 non disruptive mutants, 195 disrutpive

    piprDF = pd.DataFrame(PIPRPairs)
    piprDF = piprDF.drop_duplicates()
    print("TOTAL SIZE: ", piprDF.shape)

    # save only WT
    piprDFWT = piprDF[piprDF['type'] == 'WT']
    print("WT: ", piprDFWT.shape)
    piprDFWT = piprDFWT.drop(['type'], axis=1)
    piprDFWT.to_csv("wtFragoza.tab", sep="\t", header=True, index=False)
    # save only MT
    # save only WT
    piprDFMT = piprDF[piprDF['type'] == 'MT']
    print("MT: ", piprDFMT.shape)
    piprDFWT = piprDFMT.drop(['type'], axis=1)
    piprDFWT.to_csv("mtFragoza.tab", sep="\t", header=True, index=False)

    #number disruptive vs not
    print (piprDF.label.value_counts())
    """
    1    2073
    0     255
    """

    #mt disruptions
    print ("MT disruptions: ")
    print (piprDFMT.label.value_counts())

    """
    1    1308
    0     255
    """

    piprDF = piprDF.drop(['type'], axis=1)
    return piprDF
    #piprDF.to_csv("bothFragozaInteractions.tab", sep="\t", header=True, index=False)

#try look up interactions in orig DB
piprDF = getFragozaInteractionDFs()

print (piprDF.iloc[0])

