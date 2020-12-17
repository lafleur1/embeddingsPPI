#cleaning up DB a bit

#code to create a SQL database of PPI interactions
#making a combined DB of human PPIs w/ HuRI, HitPredict, HIPPIE, BioPlex
"""
Positives:
HIPPE: http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/ (done)
    #last updated:  Feb 14, 2019 (v2.2)
    going to want  physical association (MI:0915), direct interaction (MI:0407)
    does not have uniprot IDs/NP ids.  Has protein names. Does NOT specify isoforms interacting

HitPredict: http://www.hitpredict.org/ (done)
    Going to use physical association (MI:0915) and direct interaction only
    #last updated: 01Aug2019
    #uses uniprot IDs
HuRI:


Negatives:
Negatome 2.0 (using combined stringent): http://mips.helmholtz-muenchen.de/proj/ppi/negatome/
Non-interacting pairs (HuRI)
"""
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


#actual Sqlite code for the database now
dbName = "../interactionDB.sqlite3"
DEFAULT_PATH = os.path.join(os.path.dirname(__file__), dbName)

def db_connect(db_path=DEFAULT_PATH):
    con = sqlite3.connect(db_path)
    return con

def openMitab27File(fileName):
    #mitab27 format https://psicquic.github.io/MITAB27Format.html
    columnNames = ["Unique ID A", "Unique ID B",  "Alt ID A", "Alt ID B", "Aliases A", "Aliases B", 'Interaction Detection Methods', 'First Author', 'ID of publication', 'NCBI Tax ID A', 'NCBI Tax ID B',  'Interaction Types', 'Source Database', 'Interaction ID',  'Confidence Score', 'Complex Expansion', 'Bio role A', 'Bio role B', 'Exp role A', 'Exp role B',  'Int type A', 'Int type B', 'Xref A', 'Xref B', 'Xref Interaction',  'Annotations A', 'Annotations B', 'Annotations Int', 'NCBI Tax Host', 'Parameters Int', 'Creation Date', 'Upate Date', 'Checksum A', 'Checksum B', 'Checksum Int', 'Negative', 'Features A', 'Features B', 'Stoichiometry A', 'Stoichiometry B', 'Participant ID Method A', 'Participant ID Method B']
    huRIUnion = pd.read_csv(fileName, sep = "\t", names = columnNames )
    return huRIUnion


def justGetPairsInteractions(csvHuRI):
    csvHuRI = csvHuRI[["Unique ID A", "Unique ID B"]]
    print (csvHuRI.shape)
    #drop dups
    csvHuRI = csvHuRI.drop_duplicates()
    print (csvHuRI.shape)
    return csvHuRI

def addHuRIInteractionsPart2():
    #adds HuRI interactions to DB that have no score
    huScored = openMitab27File("./HuRI/HI-union.psi")
    print(huScored.columns)
    print("HURI TOTAL: ", huScored.shape[0])
    print (huScored)
    #huScored = huScored[huScored['Confidence Score'] == "-"] #only those not scored already, as those are already in the DB
    #print (huScored)
    conn = db_connect()
    cur = conn.cursor()
    match = 0
    updated = 0
    notMatch = 0
    duplicatePairs = []
    multIDMatches = []
    for ind, row in huScored.iterrows():
        if ind % 1000 == 0:
            print(ind, " out of ", huScored.shape[0])
            print("matches: ", match, " updates: ", updated, " not matched: ", notMatch, " number duplicated pairs: ",
                  len(duplicatePairs), " multi  matched HuRI IDs: ", len(multIDMatches))
        id_a = row["Unique ID A"].split(":")[1] #need to split off junk
        id_b = row["Unique ID B"].split(":")[1] #need to split off junk
        #just adding scored subset of HuRI for now
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HuRI", id_a))
        fetched_a = cur.fetchall()
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HuRI", id_b))
        fetched_b = cur.fetchall()
        if len(fetched_a) == 1 and len(fetched_b) == 1:
            # both are found
            # sort so that the pairs are ordered in the interactor db
            ids = []
            ids.append(fetched_a[0][3])
            ids.append(fetched_b[0][3])
            ids.sort()
            #print (ids)
            #print (score)
            # see if the combo is already in interactions
            search_interactions = "SELECT * FROM interaction WHERE protein_1=? AND protein_2=?"
            cur.execute(search_interactions, (ids[0], ids[1]))
            interaction_search = cur.fetchall()
            if len(interaction_search) == 0:
                # add new interaction
                id_map_sql = "INSERT INTO interaction (protein_1, protein_2, huri_positive) VALUES (?,?,?)"
                cur.execute(id_map_sql, (ids[0], ids[1], True))
                match += 1
                conn.commit()
            elif len(interaction_search) == 1:
                # add a HuRI score
                #updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                updateQuery = "UPDATE interaction set huri_positive=? WHERE id=?"
                cur.execute(updateQuery, (True, interaction_search[0][0]))
                updated += 1
                conn.commit()
            else:
                print("ERROR! matched more than 1 interaction......database contains duplicates....")
                duplicatePairs.append((ids[0], ids[1]))
                for matchFound in interaction_search:
                    # add HuRI score to each duplicate
                    updateQuery = "UPDATE interaction set huri_positive=? WHERE id=?"
                    cur.execute(updateQuery, (True, matchFound[0]))
                    updated += 1
                    conn.commit()
                    #updated += 1
        else:
            print("FOUND: ", len(fetched_a), " protein ids interactor a ", len(fetched_b), " protein ids interactor b")
            notMatch += 1
            if len(fetched_a) > 1:
                print("matches mult. ids....")
                multIDMatches.append(id_a)
            elif len(fetched_b) > 1:
                print("matches mult ids......")
                multIDMatches.append(id_b)
    print("total: ", huScored.shape)
    print("matches: ", match)
    print(" not matches: ", notMatch)
    print ("number bad ids....", len(multIDMatches))
    print("____________________________________________________")
    for id in multIDMatches:
        print (id)
    print ("____________________________________________________")
    return duplicatePairs

def lookAtPairs(pairs):
    #making sure HuRI is properly inserted into DB again
    pairsRetrieved = []
    conn = db_connect()
    cur = conn.cursor()
    foundPairs = 0
    notFound = 0
    intNotPres = 0
    intPres = 0
    intPresMult = 0
    for ind, row in pairs.iterrows():
        id_a = row["Unique ID A"].split(":")[1]  # need to split off junk
        id_b = row["Unique ID B"].split(":")[1]  # need to split off junk
        # just adding scored subset of HuRI for now
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HuRI", id_a))
        fetched_a = cur.fetchall()
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HuRI", id_b))
        fetched_b = cur.fetchall()
        if len(fetched_a) == 1 and len(fetched_b) == 1:
            #print (fetched_a, fetched_b)
            foundPairs += 1
            #prep to add to interactions
            ids = []
            ids.append(fetched_a[0][3])
            ids.append(fetched_b[0][3])
            ids.sort()
            #try to see if it already exists
            search_interactions = "SELECT * FROM interaction WHERE protein_1=? AND protein_2=?"
            cur.execute(search_interactions, (ids[0], ids[1]))
            interaction_search = cur.fetchall()
            if len(interaction_search) == 0:
                intNotPres += 1
            elif len(interaction_search) == 1:
                intPres += 1
            elif len(interaction_search) > 1:
                intPresMult += 1
        else:
            print ("ERROR!")
            notFound += 1
            print (row)

    print (foundPairs, " ", notFound, " ", foundPairs + notFound)
    print ("those found: ", intNotPres, " ", intPres, " ", intNotPres)
    """
    appears to be properly handed in DB 
    96099   1789   97888
    those found:  201   95212   201
    """

if __name__ == '__main__':
    mitab = openMitab27File("./HuRI/HI-union.psi")#addHuRIInteractionsPart2()
    #print (mitab.iloc[0])
    combos = justGetPairsInteractions(mitab)
    #addUniqueInteractions(combos)
    #mitab.to_csv("demo.csv")
    lookAtPairs(combos)