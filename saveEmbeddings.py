#save embeddings for each fasta in the folder

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

#bert_model =  ProteinBertModel.from_pretrained('bert-base')
bert_model = None
bert_tokenizer = None
#bert_tokenizer = TAPETokenizer(vocab='iupac') 
unirep_model = UniRepModel.from_pretrained('babbler-1900')
unirep_tokenizer = TAPETokenizer(vocab='unirep')


def db_connect(db_path=DEFAULT_PATH):
    con = sqlite3.connect(db_path)
    return con

def getTransformerMeanSequenceOutput(seq, model = bert_model, tokenizer = bert_tokenizer):
    #from github: pooled_output is *not* trained for the transformer
    #model = ProteinBertModel.from_pretrained('bert-base')
    #tokenizer = TAPETokenizer(vocab='iupac')  # iupac is the vocab for TAPE models, use unirep for the UniRep model
    token_ids = torch.tensor([tokenizer.encode(seq)])
    output = model(token_ids)
    #average embedding for the bert model, average along sequence as recommendd
    #output shape is (1,768)
    return torch.mean(output[0],1)

def getUniRepPooledSequence(sequence, model = unirep_model, tokenizer = unirep_tokenizer):
    #return
    #model = UniRepModel.from_pretrained('babbler-1900')
    #tokenizer = TAPETokenizer(vocab='unirep')
    token_ids = torch.tensor([tokenizer.encode(sequence)])
    output = model(token_ids)
    return output[1]


def getCDHitToEmbed(clusterFile, outputDir):
    output = []
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
            seqs.append(fetched[0][1])
            output.append((fetched[0][0], fetched[0][1]))
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print ("found: ", match)
    print ("not found: ", notMatch)
    return seqs, output 

def createReps(clusterFile, outputDir):
    seqs, output = getCDHitToEmbed(clusterFile, outputDir)
    '''
    for i in range(0, len(output)):# in output:
        pair = output[i]
        print ("ON: ", i, " OUT OF ", len(output))
        print("len seq: ", len(pair[1]))
        meanS1 = getTransformerMeanSequenceOutput(pair[1]).detach().numpy()
        #print(meanS1.shape)
        #print (type(meanS1))
        #save npy file to the folder 
        np.save("./bert_base_embeddings/" + str(pair[0]) + ".npy", meanS1)
        print ("bert done")
    '''
    output = output[::-1]
    for i in range(0, len(output)):# in output:
        pair = output[i]
        print ("ON: ", i, " OUT OF ", len(output))
        print("len seq: ", len(pair[1]))
        if not path.exists("./unirep_embeddings/" + str(pair[0]) + ".npy"):
            dummyFile = np.array([0])
            np.save("./unirep_embeddings/" + str(pair[0]) + ".npy", dummyFile) #save temp blank file here to run mult copies of the script at once 
            uniRepPool = getUniRepPooledSequence(pair[1]).detach().numpy()
            #print(uniRepPool.shape)
            #print (type(uniRepPool))
            np.save("./unirep_embeddings/" + str(pair[0]) + ".npy", uniRepPool)
            print ("unirep done")
        else:
            print ("exists")
            
        
createReps("embedPPI_50_2000_70_Cutoff.clstr", "./embedSeqs/")
