#setting up pytorch & trying the demos
import torch
from tape import ProteinBertModel, TAPETokenizer
#unirep
from tape import UniRepModel
#getting the sweet, sweet < 2K AA proteins to use for the TAPE embeddings
import string
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyensembl
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import re
from urllib.request import urlopen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from io import TextIOWrapper
import time
import os.path
import matplotlib.pyplot as plt
from os import path
from queue import PriorityQueue
from itertools import combinations
from numpy.random import choice
import random
import sqlite3
#from venn import venn
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#actual Sqlite code for the database now
dbName = "interactionDB.sqlite3"
DEFAULT_PATH = os.path.join(os.path.dirname(__file__), dbName)

def db_connect(db_path=DEFAULT_PATH):
    con = sqlite3.connect(db_path)
    return con

def getSequences():
    #return the contents of the sequence table
    conn = db_connect()
    cur = conn.cursor()
    with conn:
        cur.execute("SELECT * FROM proteins")
        return cur.fetchall()

def getAllInteractions():
    # return the contents of the sequence table
    conn = db_connect()
    cur = conn.cursor()
    with conn:
        cur.execute("SELECT * FROM interaction")
        return cur.fetchall()



def getSubsetInteractions(idList):
    # return all interactions where both partners in idList
    conn = db_connect()
    cur = conn.cursor()
    cur.execute("SELECT * FROM interaction")
    ints = cur.fetchall()
    subset = []
    for i in range(0, len(ints)):
        int_potential = ints[i]
        if i%1000 == 0:
            print ("ON: ", i, " OUT OF ", len(ints))
            print (len(subset))
        #check each interaction partners if they are in the list of ids
        #print (int_potential)
        if int_potential[1] in idList and int_potential[2] in idList:
            subset.append(int_potential)
    return subset

def resetAllSequenceLens():
    #getting all sequences <= 2000 AA long
    allseqs = getSequences()
    ids = [x[0] for x in allseqs]
    actualSeqs = [x[1] for x in allseqs]
    lens = [len(x) for x in actualSeqs]
    conn = db_connect()
    cur = conn.cursor()
    #set all sequence lengths
    for i in range(0, len(ids)):
        updateQuery = "UPDATE proteins SET length=? WHERE id=?"
        cur.execute(updateQuery, (lens[i], ids[i]))
        conn.commit()

def proteinsInDBtoFastaForCDHit(lower, upper, saveName):
    #convertes the SQL table to a fasta with my_ids as the name per sequence to run 40% CD-HIT on
    #distinctProQuery = "SELECT DISTINCT sequence FROM proteins"
    conn = db_connect()
    cur = conn.cursor()
    #query = cur.execute(distinctProQuery)
    #fetched = query.fetchall()
    #confirm that total size of db is same as unique number of seqs
    allQuery = "SELECT * FROM proteins WHERE length BETWEEN " + str(lower) + " AND " + str(upper) #trRosetta trained on sequences <= upper and >= lower
    allQuery = cur.execute(allQuery)
    allFetched = allQuery.fetchall()
    #if len(fetched) != len(allFetched):
     #   print ("stopping, duplicate sequences in teh table")
    #else:
    print ("writing all sequences to a fasta")
    print (len(allFetched))
    #make into a large fasta file w/ CD-HIT
    with open(saveName, "w") as f:
        for result in allFetched:
            f.write(">" + str(result[0]) + "\n" + result[1] + "\n" )


def getClusterRepIDs(clusterFile):
    #each rep has a * by the name
    reps = []
    with open(clusterFile, "r") as f:
        lines = [x.strip() for x in f.readlines()]
        for l in lines:
            if "*" in l:
                half = l.split(">")[1].translate(str.maketrans('', '', string.punctuation))
                half = half.replace(" ","")
                reps.append(half)
    #change to ints
    reps = [int(x) for x in reps]
    return reps

#todo: extract ids of CD-HIT clusters
def getCDHIT(clusterFile, outputDir):
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    #each rep has a * by the name
    reps = []
    with open(clusterFile, "r") as f:
        lines = [x.strip() for x in f.readlines()]
        for l in lines:
            if "*" in l:
                half = l.split(">")[1].translate(str.maketrans('', '', string.punctuation))
                half = half.replace(" ","")
                reps.append(half)
    #get sequences from db
    conn = db_connect()
    cur = conn.cursor()
    match = 0
    notMatch = 0
    seqs = []
    print ('total len reps: ', len(reps))
    for id in reps:
        cur.execute("SELECT * FROM proteins WHERE id=?", (id,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            # add id to the my_ids table
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            seqs.append(fetched[0][1])
            match += 1
            with open(outputDir + id + ".fasta", "w") as f:
                f.write(">" + str(id) + "\n" + fetched[0][1] + "\n")
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print ("found: ", match)
    print ("not found: ", notMatch)
    return seqs

def getTransformerMeanSequenceOutput(seq):
    #from github: pooled_output is *not* trained for the transformer
    model = ProteinBertModel.from_pretrained('bert-base')
    tokenizer = TAPETokenizer(vocab='iupac')  # iupac is the vocab for TAPE models, use unirep for the UniRep model
    token_ids = torch.tensor([tokenizer.encode(seq)])
    output = model(token_ids)
    #average embedding for the bert model, average along sequence as recommendd
    #output shape is (1,768)
    return torch.mean(output[0],1)

def getUniRepPooledSequence(sequence):
    #return
    model = UniRepModel.from_pretrained('babbler-1900')
    token_ids = torch.tensor([tokenizer.encode(sequence)])
    output = model(token_ids)
    return output[1]


def retrieveOrigNames(namesDone):
    #retrieve count of proteins in HuRI
    con = db_connect()
    cur = con.cursor()
    huRIMappedNames = []
    for i in range(0, len(namesDone)):
        name = namesDone[i]
        if i % 1000 == 0:
            print ("ON ", i , " OUT OF ", len(namesDone))
        lookupQ = cur.execute("SELECT * FROM id_map WHERE database=? AND my_id=?", ("HuRI", name))
        fetched = lookupQ.fetchall()
        if len(fetched) != 0:
            #print ("found match in HuRI!!!")
            huRIMappedNames.append(fetched[0][2])
    #print (len(huRIMappedNames))
    return huRIMappedNames

def overlapTwoLists(list1, list2):
    intersection = set(list1).intersection(set(list2))
    return list(intersection)

def lookAtPositivesAndNegativesCurrent(listIDs):
    # looking at the number of positive examples with the cluster reps
    cdHitSeqs = getClusterRepIDs("embedPPI_50_2000_70_Cutoff.clstr")  # sequence IDs in the DB
    # get all interactions
    subsetGoodFor = getSubsetInteractions(cdHitSeqs)
    positives = [x for x in subsetGoodFor if x[8] != 1]
    print(len(positives))  # 279775
    negatives = [x for x in subsetGoodFor if x[8] == 1]
    print(len(negatives))  # 918
    return positives, negatives


def lookAtIntersection(cdHitSeqs):
    #
    print ('start ret')
    huRIMappedNames = retrieveOrigNames(cdHitSeqs)
    print(len(huRIMappedNames))
    in05 = open("./ppiDB/HuRI/HuRI_05_proteins.txt", "r")
    in05 = [x.strip() for x in in05.readlines()]
    intersect = overlapTwoLists(huRIMappedNames, in05)
    print(len(intersect))
    print(len(in05))
    in14 = open("./ppiDB/HuRI/HuRI_14_proteins.txt", "r")
    in14 = [x.strip() for x in in14.readlines()]
    intersect = overlapTwoLists(huRIMappedNames, in14)
    print(len(intersect))
    print(len(in14))
    in19 = open("./ppiDB/HuRI/HuRI_HuRI_proteins.txt", "r")
    in19 = [x.strip() for x in in19.readlines()]
    intersect = overlapTwoLists(huRIMappedNames, in19)
    print(len(intersect))
    print(len(in19))
    # retrieve HuRI associated names to construct HuRI negatives
    # print (namesDone)


#representative sequence clustering:
#1) extract all appropriate length sequences to a fasta file
#proteinsInDBtoFastaForCDHit(50, 2000, "embedPPI_50_2000.fasta")
#2) Run CD-HIT and get actual representative sequences post clustering
#then run CD-HIT with: cd-hit -i embedPPI_50_2000.fasta -o embedPPI_50_2000_70_Cutoff -c 0.7 to generate cluster
#3) Extract cluster representatives (16614 proteins between 50 & 20000 AA at 70% identity threshold)
# 4) Save each as npy file using the embedding methods


#looking at the number of positive examples with the cluster reps
#cdHitSeqs = getClusterRepIDs("embedPPI_50_2000_70_Cutoff.clstr") #sequence IDs in the DB
#lookAtIntersection(cdHitSeqs)

