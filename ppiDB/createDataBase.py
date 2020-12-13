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

def isHighConfidenceHitPredict( hitPredictInteractionScore):
    #returns if is high confidence or not based on HitPredict score cutoffs
    #If: annotation score >= 0.5 or method score >= 0.485 or interaction score > 0.28 is considered high confidence
    if hitPredictInteractionScore > 0.28:
        return True
    else:
        return False

def isHighConfidenceHIPPIE(hippieScore):
    #returns if the interaction is high confidence based on HIPPIE cutoffs
    if hippieScore >= 0.73:
        return True
    else:
        return False


def openHitPredict():
    #opens the hit predict mitab adn returns as a pandas DF
    """
    Index(['Interactor A ID', 'Interactor B ID', 'Interactor A Name',
       'Interactor B Name', 'Interactor A Alias', 'Interactor B Alias',
       'Interaction detection methods', 'First author', 'Publications',
       'Interactor A Taxonomy', 'Interactor B Taxonomy', 'Interaction type',
       'Source database', 'Interaction identifier', 'Confidence score'],
      dtype='object')
    """
    hitPredict =  pd.read_csv("./HitPredict/HitPredict_H_sapiens_interactions_MITAB-2.5", skiprows = 5, sep= "\t")
    #check all are 9606 (both partners are human)
    hitPredict = hitPredict[(hitPredict['Interactor A Taxonomy'] == 9606) &(hitPredict['Interactor B Taxonomy'] == 9606) ].copy().reset_index(drop = True)
    #remove uniprotkb: from start of all ids
    hitPredict['UniprotID_A'] = hitPredict.apply(lambda row: row['Interactor A ID'].split(":")[1], axis = 1)
    hitPredict['UniprotID_B'] = hitPredict.apply(lambda row: row['Interactor B ID'].split(":")[1], axis=1)
    # keep only those with MI:0915 (physical association) and MI:0407 (direction interaction)
    #print ("HitPredict all: ", hitPredict.shape)
    #hitPredict = hitPredict[(hitPredict['Interaction type'].str.contains("psi-mi:MI:0915")) | (hitPredict['Interaction type'].str.contains("psi-mi:MI:0407"))].copy().reset_index(drop = True)
    # drop pubmed stuff, alias columns (not used)
    #print ("HitPredict direction assoc or physical assoc only: ", hitPredict.shape)
    hitPredict = hitPredict.drop(
        ['Interactor A Alias', 'Interactor B Alias', 'Publications', 'First author', 'Interactor A Taxonomy',
         'Interactor B Taxonomy','Interactor A ID' , 'Interactor B ID'], axis=1)
    obsoleteIDs = ["D3DNK7", "A0A087WXL1", "B5MC15", "A0A087WUM8", "Q5JXX5", "A6NE70", "Q5SQP8", 'A2A2G4', 'J3KNN5', 'A0A0A0MTS2', 'A0A024RDG1', 'A0A024R419', 'J3QT77', 'A0A024R912', 'A0A024R473', 'A0A024R2M7', 'H0YKG7', 'A0A024R2T6', 'A0A024RD93', 'A0A024R4H0', 'B1AKJ5', 'A0A024R3Z6', 'A0A024R497', 'K7EK91']
    #remove all interactions which are obsolete in UniProt as of 7/6/20
    print ("Remove all interactions which are obsolete in UniProt as of 7/6/20")
    print ("pre-removal: ", hitPredict.shape)
    hitPredict = hitPredict[~(hitPredict['UniprotID_A'].isin(obsoleteIDs)) & ~(hitPredict['UniprotID_B'].isin(obsoleteIDs)) ]
    print ("post-removal: ", hitPredict.shape)
    hitPredict['High Confidence'] = hitPredict.apply(lambda row: isHighConfidenceHitPredict(row['Confidence score']), axis = 1)
    return hitPredict

def getAllHitPredictUniprotKBIDs():
    #save all uniprot IDs from HitPredict to get sequences for
    hitPredict = openHitPredict()
    geneA = hitPredict['UniprotID_A'].to_list()
    geneB = hitPredict['UniprotID_B'].to_list()
    allGenes = set(geneA + geneB)
    print ("total proteins: ", len(allGenes))
    file = open("allHitPredictUniprotIds.txt", "w")
    for name in allGenes:
        file.write(name + "\n")
    file.close()
    return allGenes

def openHIPPIE():
    #open HIPPIE ppis
    """
    Index(['ID Interactor A', 'ID Interactor B', 'Alt IDs Interactor A',
       'Alt IDs Interactor B', 'Aliases Interactor A', 'Aliases Interactor B',
       'Interaction Detection Methods', 'Publication 1st Author',
       'Publication Identifiers', 'Taxid Interactor A', 'Taxid Interactor B',
       'Interaction Types', 'Source Databases', 'Interaction Identifiers',
       'Confidence Value', 'Presence In Other Species',
       'Gene Name Interactor A', 'Gene Name Interactor B'],
      dtype='object')
    """
    hippie = pd.read_csv("./HIPPIE/HIPPIE-current.mitab.txt", sep = "\t")
    #filter for only human
    #print("hippie size: ", hippie.shape)
    hippie = hippie[(hippie['Taxid Interactor A'] == 'taxid:9606(Homo sapiens)') & (hippie['Taxid Interactor B'] == 'taxid:9606(Homo sapiens)')].copy().reset_index(drop = True)
    #print ("hippie size only human: ", hippie.shape)
    #is in vivo (0493) or MI:0018(two hybrid)
    #fill in high confidence boolean
    hippie['High Confidence'] = hippie.apply(lambda row: isHighConfidenceHIPPIE(row['Confidence Value']), axis = 1)
    #print ("only high confidence: " )
    #print (hippie['High Confidence'].value_counts())
    #hippie = hippie[(hippie['Interaction Detection Methods'].str.contains("two hybrid")) | (hippie['Interaction Detection Methods'].str.contains("two hybrid"))].copy().reset_index(drop = True)
    #print ("only those with two hybrid as a experimental detection method: ", hippie.shape)
    # remove columns
    #print(hippie['High Confidence'].value_counts())
    #drop those not associated with a gene
    hippie = hippie[(hippie['ID Interactor A'].str.contains("entrez gene")) & (hippie['ID Interactor B'].str.contains("entrez gene"))].copy().reset_index(drop=True)
    hippie['Gene ID A'] = hippie['ID Interactor A'].apply(lambda x: x.split(":")[1])
    hippie['Gene ID B'] = hippie['ID Interactor B'].apply(lambda x: x.split(":")[1])
    #getting uniprot kb protein names
    hippie = hippie[(hippie['Alt IDs Interactor A'].str.contains("uniprotkb")) & (
        hippie['Alt IDs Interactor B'].str.contains("uniprotkb"))].copy().reset_index(drop=True)
    hippie['Uniprot Name A'] = hippie['Alt IDs Interactor A'].apply(lambda x: x.split(":")[1])
    hippie['Uniprot Name B'] = hippie['Alt IDs Interactor B'].apply(lambda x: x.split(":")[1])
    #print (hippie.shape)
    hippie = hippie.drop(
        ['Taxid Interactor A', 'Taxid Interactor B', 'Publication Identifiers', 'Presence In Other Species',
         'Interaction Types', 'Interaction Detection Methods','ID Interactor A','ID Interactor B', 'Alt IDs Interactor A', 'Alt IDs Interactor B', 'Aliases Interactor A', 'Aliases Interactor B', 'Interaction Identifiers', 'Publication 1st Author', 'Source Databases'], axis=1)
    return hippie

def getAllHIPPIEUniprotNames():
    hippie = openHIPPIE()
    genesA = hippie['Gene ID A'].to_list()
    genesB = hippie['Gene ID B'].to_list()
    allGenes = set(genesA + genesB)
    namesA = hippie['Uniprot Name A'].to_list()
    namesB = hippie['Uniprot Name B'].to_list()
    names = set(namesA + namesB)
    print ("total names: ", len(names))
    file = open("hippieUniprotNames.txt", "w")
    for name in names:
        file.write(name + "\n")
    file.close()
    return allGenes, names

def extractRelevantInfoHitPredict():
    #creating information for mapping names for making SQL database
    """
    Index(['Entry', 'HitPredictID', 'Cross-reference (GeneID)', 'Entry name',
       'Protein names', 'Length', 'Status', 'Sequence'],
      dtype='object')
    """
    mappings = [] #tuples of uniprot id, hitpredict id pairs
    uniprotHitPredict = pd.read_csv("./HitPredict/uniprotHitPredict.tab", sep = "\t")
    doubleMaps = uniprotHitPredict[uniprotHitPredict['HitPredictID'].str.contains(",")].copy().reset_index(drop = True)
    singleMaps = uniprotHitPredict[~uniprotHitPredict['HitPredictID'].str.contains(",")].copy().reset_index(drop = True)
    print ("double mapped: ", doubleMaps.shape) #9 entries
    for ind, row in doubleMaps.iterrows():
        allHitPredictIDs = row['HitPredictID'].split(",")
        for id in allHitPredictIDs:
                mappings.append((row['Entry'], id))
    print (mappings)
    print ("single mapped: ", singleMaps.shape) #20498
    singleMaps['mapping'] = list(zip(singleMaps['Entry'], singleMaps['HitPredictID']))
    mappings = list(set(mappings + singleMaps['mapping'].to_list()))
    print (len(mappings)) #total is 20559, 100% mappings recovered
    #justSeqInfoHitPredict = uniprotHitPredict.drop(['HitPredictID'], axis = 1)
    return mappings, uniprotHitPredict


def extractRelevantInfoHIPPIE():
    #extract mappings from HIPPIE protein names to uniprot entry ids, throw away non-essential information
    """
    Index(['HIPPIENames', 'Entry',
       'yourlist:M20200706A94466D2655679D1FD8953E075198DA8141F74D',
       'Cross-reference (GeneID)', 'Entry name', 'Protein names', 'Length',
       'Status', 'Sequence'],
      dtype='object')
    """
    #manually removing names which are duplicated unecessarily in uniprot HIPPIE retrieval
    badNames = open("./HIPPIE/badNamesHIPPIE.txt", "r")
    badUniprotIDs = [x.strip() for x in badNames.readlines()]
    badNames.close()
    #set up mappings dictionary
    mappings = []
    mappingsDict = {'hippieName':[], 'entry':[], 'geneID':[]}
    uniprotHIPPIE = pd.read_csv("./HIPPIE/uniprotHIPPIEV2.tab", sep="\t") #17229 rows
    #drop duplicate sequences
    uniprotHIPPIE = uniprotHIPPIE.drop_duplicates(subset = "Sequence")
    #remove those which are in the bad duplicates
    print (uniprotHIPPIE.shape)
    uniprotHIPPIE = uniprotHIPPIE[~uniprotHIPPIE['Entry'].isin(badUniprotIDs)].copy().reset_index(drop = True)
    print (uniprotHIPPIE.shape)
    print ("HIPPIE uniprot: ", uniprotHIPPIE.shape)
    doubleMaps = uniprotHIPPIE[uniprotHIPPIE['HIPPIENames'].str.contains(",")].copy().reset_index(drop=True)
    doubleMaps.to_csv("./HIPPIE/doubleMapHIPPIE.csv")
    print ("multi mapped: ", doubleMaps.shape)
    for ind, row in doubleMaps.iterrows():
        allNames = row['HIPPIENames'].split(",")
        #print (len(allNames))
        for n in allNames:
            mappings.append((row['Entry'], n))
            mappingsDict['hippieName'].append(n)
            mappingsDict['entry'].append(row['Entry'])
            mappingsDict['geneID'].append(row['Cross-reference (GeneID)'])
    singleMaps = uniprotHIPPIE[~uniprotHIPPIE['HIPPIENames'].str.contains(",")].copy().reset_index(drop=True)
    print ("single mapped: ", singleMaps.shape)
    singleMaps['mapping'] = list(zip(singleMaps['Entry'], singleMaps['HIPPIENames']))
    mappingsDict['hippieName'] = mappingsDict['hippieName'] + singleMaps['HIPPIENames'].to_list()
    mappingsDict['entry'] = mappingsDict['entry'] + singleMaps['Entry'].to_list()
    mappingsDict['geneID'] = mappingsDict['geneID'] + singleMaps['Cross-reference (GeneID)'].to_list()
    mappings = list(set(mappings + singleMaps['mapping'].to_list()))
    #print(len(mappings))  # total is 20559, 100% mappings recovered
    getAllHIPPIEUniprotNames()
    #print (len(mappingsDict))
    toDF = pd.DataFrame(mappingsDict)
    names = toDF['hippieName'].value_counts().to_frame()
    names.to_csv("hippieNames.csv")
    toDF.to_csv("hippieMappings.csv")
    #keep only sequence information for the entries which match with a HIPPIE protein
    entriesListMapped = toDF['entry'].to_list()
    uniprotHIPPIE = uniprotHIPPIE[uniprotHIPPIE['Entry'].isin(entriesListMapped)].copy().reset_index(drop = True)
    #drop mylists
    #uniprotHIPPIE = uniprotHIPPIE.drop(["HIPPIENames", 'yourlist:M20200707216DA2B77BFBD2E6699CA9B6D1C41EB2968714N'], axis = 1)
    print ("remaining sequences: ", uniprotHIPPIE.shape)
    return toDF, uniprotHIPPIE

def manuallySelectCorrectUniprotsForDoubledHIPPIENames():
    extractRelevantInfoHIPPIE()
    print("DONE")
    #open hippie mappings and fix column names
    hippie = openHIPPIE()
    # print (hippie.iloc[0])
    mappingsDict = pd.read_csv("hippieMappings.csv", sep=",")
    justNamesBad = pd.read_csv("hippieNames.csv", sep=",")
    justNamesBad = justNamesBad[justNamesBad['matches'] > 1]
    print(justNamesBad)
    for ind, row in justNamesBad.iterrows():
        print("--------------------------------------------------------dddddddddd")
        print(row)
        print(mappingsDict[mappingsDict['hippieName'] == row["hippieName"]])
        print(hippie[(hippie['Uniprot Name A'] == row["hippieName"]) | (hippie['Uniprot Name B'] == row["hippieName"])])
        input()

def openMitab27File(fileName):
    #mitab27 format https://psicquic.github.io/MITAB27Format.html
    columnNames = ["Unique ID A", "Unique ID B",  "Alt ID A", "Alt ID B", "Aliases A", "Aliases B", 'Interaction Detection Methods', 'First Author', 'ID of publication', 'NCBI Tax ID A', 'NCBI Tax ID B',  'Interaction Types', 'Source Database', 'Interaction ID',  'Confidence Score', 'Complex Expansion', 'Bio role A', 'Bio role B', 'Exp role A', 'Exp role B',  'Int type A', 'Int type B', 'Xref A', 'Xref B', 'Xref Interaction',  'Annotations A', 'Annotations B', 'Annotations Int', 'NCBI Tax Host', 'Parameters Int', 'Creation Date', 'Upate Date', 'Checksum A', 'Checksum B', 'Checksum Int', 'Negative', 'Features A', 'Features B', 'Stoichiometry A', 'Stoichiometry B', 'Participant ID Method A', 'Participant ID Method B']
    huRIUnion = pd.read_csv(fileName, sep = "\t", names = columnNames )
    return huRIUnion

def getAllHuRIProteinNames(hu, uniprotTxtFileName  = "HuRIUniprotNames.txt", ensemblTxtFileName = "EnsemblUniprotNames.txt"):
    #extracts all uniprot and ensembl ids for proteins
    #input of huri file to get ids from
    #saves ids to seprate text files for uniprot Ids and ensembl ids
    #hu = openMitab27File(fileName)
    names = {'Uniprot':[], 'Ensembl':[]}
    allNames = list(set(hu['Unique ID A'].to_list() + hu['Unique ID B'].to_list()))
    #manually separate
    for name in allNames:
        if "uniprotkb" in name:
            names['Uniprot'].append(name.split(":")[1])
        elif 'ensembl' in name:
            names['Ensembl'].append(name.split(":")[1])
        else:
            print ("ERROR!  Unknown name...", name)
    print ("total uniprot names: ", len(names['Uniprot']))
    print ("total ensembl names: ", len(names['Ensembl']))
    upFile = open(uniprotTxtFileName, "w")
    for name in names['Uniprot']:
        upFile.write(name + "\n")
    upFile.close()
    enFile = open(ensemblTxtFileName, "w")
    for name in names["Ensembl"]:
        enFile.write(name + "\n")
    enFile.close()

def extractOnlyScoredPPIsHuRI2019():
    hu2019 = openMitab27File("./HuRI/HuRI.psi") #2019 Luck HuRI dataset has most PPIs associated with a confidence score
    scoredHu = hu2019[~(hu2019['Confidence Score'] == "-")].copy().reset_index(drop=True)
    return scoredHu


def manualEnsemblViewing():
    #manually sorted the ensembl ids, as the ensembl-uniprot database links are not udpated and id mapping tool does not work
    ensemblNames = pd.read_csv("ensemblToLocation.csv")
    inFolder = ensemblNames[ensemblNames['Uniprot'] == "FOLDER"]
    print("In folder: ", inFolder.shape)  # folder 18
    deleted = ensemblNames[ensemblNames['Uniprot'] == "DELETED"]
    print("Deleted: ", deleted.shape)  # 6 deleted
    uniprot = ensemblNames[~(ensemblNames['Uniprot'] == "FOLDER") & ~(ensemblNames['Uniprot'] == "DELETED")]
    print("uniprot ids: ", uniprot.shape)  # 32 uniprot ids
    uniprotNamesEnsembl = uniprot['Uniprot'].to_list()
    openF = open("HuRIEnsemblUniprot.txt", "w")
    openF.write("\n".join(uniprotNamesEnsembl))
    openF.close()

def negatomeIDs():
    negatome = pd.read_csv("Negatome2_combinedStringent.txt", sep = "\t", header= None)
    names = list(set(negatome[0].to_list() + negatome[1].to_list()))
    print ("total: ", len(names))
    openF = open("negatomeUniprotIDs.txt", "w")
    openF.write("\n".join(names))
    openF.close()


#actual Sqlite code for the database now
dbName = "interactionDB.sqlite3"
DEFAULT_PATH = os.path.join(os.path.dirname(__file__), dbName)

def db_connect(db_path=DEFAULT_PATH):
    con = sqlite3.connect(db_path)
    return con

def setupDB():
    # sequence structure is going to be id (int) sequence (string)
    protein_sequence = """
        CREATE TABLE proteins (
            id integer PRIMARY KEY,
            sequence text NOT NULL ) """
    #idMap is going to be id in database, database source, and foreign key in protein table
    id_map = """
        CREATE TABLE id_map (
        id integer PRIMARY KEY,
        database text NOT NULL,
        database_name text NOT NULL,
        my_id integer,
        FOREIGN KEY (my_id) references protein_sequence (id),
        UNIQUE(database,database_name))
    """
    #interaction is going to be my id 1, my id 2 (ordered by increasing number so that 1->2 and 2->1 is not an issue)
    #boolean if is positive or is negative
    #HIPPIE_score
    #HitPredict score
    #HuRI score
    interactions = """
            CREATE TABLE interaction (
                id integer PRIMARY KEY,
                protein_1 int,
                protein_2 int,
                hippie_score real,
                hippie_high_confidence boolean,
                hitpredict_score real,
                hitpredict_high_confidence boolean,
                huri_score real,
                negatome_interaction boolean,
                FOREIGN KEY (protein_1) references protein_sequence (id),
                FOREIGN KEY (protein_2) references protein_sequence (id))
            """
    #create connection and tables
    con = db_connect()
    cur = con.cursor()
    cur.execute(protein_sequence)
    cur.execute(id_map)
    cur.execute(interactions)
    con.commit()

def getListAllAlternativeIsoforms():
    #gets HuRI isoforms needed from the Uniprot isoform fasta
    allNames = open("./HuRI/HuRIUniprotNames.txt", 'r')
    allNames = [x.strip() for x in allNames.readlines()]
    pullTable = []
    pullIsoforms = []
    for name in allNames:
        #if no -# or -1 put in pullTable, uses cannonical isoform for the id
        if "-" in name:
            isoform = name.split("-")
            if isoform[1] == "1":
                pullTable.append(name)
            else:
                pullIsoforms.append(name)
        else:
            pullTable.append(name)
    print ("cannonical isoforms: ", len(pullTable))
    #pull these from the table
    print ("noncannonical isoforms: ", len(pullIsoforms))
    return pullTable, pullIsoforms

def openCannonicalAndIsoformFastas():
    #takes a list of isoform ids and pulls the sequence from teh Uniprot alternative isoforms fasta file
    #turns into dictionary with key:sequence entries
    cannonical, noncannonical = getListAllAlternativeIsoforms()
    keeperSeqs = []
    foundIds = [] #mappning of tuples from HuRI to uniprot ids
    idListIntermed  = []
    #print (len(noncannonical))
    with open("./HuRI/cannonicalAndIsoformsHuRI.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            if actualID in noncannonical:
                keeperSeqs.append(str(record.seq))
                foundIds.append((actualID, str(record.seq)))
                idListIntermed.append(actualID)
            elif actualID in cannonical:
                keeperSeqs.append(str(record.seq))
                foundIds.append((actualID, str(record.seq)))
                idListIntermed.append(actualID)
    #print ("found: ", len(keeperSeqs))
    #figure out what is not there
    foundSet = set(idListIntermed)
    notFound = list(set(noncannonical).difference(foundSet))
    notFoundCannon = list(set(cannonical).difference(foundSet))
    print ("len not found cannonical: ", len(notFoundCannon))
    print("len not found noncannonical: ", len(notFound))
    justOnesNotFoundCannon = [x for x in notFoundCannon if "-1" in x]
    print ("len not found cannonical w/ -1: ", len(justOnesNotFoundCannon))
    #print (len(notFound))
    with open("./HuRI/isoformFastas.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            #print ("id: ", actualID)
            for missing in notFound:
                stub = missing.split("-")[0]
                if stub == actualID:
                    keeperSeqs.append(str(record.seq))
                    foundIds.append((missing, str(record.seq)))
                    idListIntermed.append(missing)
    #print("found: ", len(keeperSeqs))
    #print ("found: ", len(set(keeperSeqs)))
    #print ("ids: ", len(foundIds))
    #getting -1 isoforms not found
    with open("./HuRI/cannonicalAndIsoformsHuRI.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            for missing in justOnesNotFoundCannon:
                stub = missing.split("-")[0]
                if stub == actualID:
                    #looking for -1 isoforms left over which are not
                    keeperSeqs.append(str(record.seq))
                    foundIds.append((missing, str(record.seq)))
                    idListIntermed.append(missing)
    foundSet = set(foundIds)
    notFound = list(set(noncannonical).difference(foundSet))
    print ("total ids: ", len(cannonical) + len(noncannonical))
    print ("Found: ", len(foundIds))
    print ("found seqs: ", len(keeperSeqs))
    totalIds = set(cannonical + noncannonical)
    notFound = list(totalIds.difference(set(idListIntermed)))
    print ("total not found: ", len(notFound))
    #recovered 25 more from Unipro, some ids were compeltely changed
    huRI = pd.read_csv("./HuRI/replaceUniprotNames.tab", sep = "\t")
    keepers = huRI[huRI['yourlist'].isin(notFound)]
    #strip out both lists and add
    keepers['mappings'] = list(zip(keepers['yourlist'], keepers['Sequence']))
    foundIds = foundIds + keepers['mappings'].to_list()
    keeperSeqs = keeperSeqs + keepers['Sequence'].to_list()
    idListIntermed = idListIntermed + keepers['yourlist'].to_list()
    print("total ids: ", len(cannonical) + len(noncannonical))
    print("Found: ", len(set(idListIntermed)))
    print("found seqs: ", len(set(keeperSeqs)))
    return foundIds, keeperSeqs

def combinedSequences():
    #open and combine by unique sequences for HIPPIE, HitPredict, HuRI (Scored subset), and Negatome which could be recovered
    #HitPredict
    mappings, hitPredict  = extractRelevantInfoHitPredict()
    print (mappings)
    #hitPredict['Sequence'].to_list() #list of unique proteins HitPredict
    #HIPPIE
    mappingsDF, hippie = extractRelevantInfoHIPPIE()
    #hippie['Sequence'].to_list()
    #HuRI
    huriIds, keeperSeqs = openCannonicalAndIsoformFastas()
    allSequences = list(set(hitPredict['Sequence'].to_list() + hippie['Sequence'].to_list() + keeperSeqs))
    print (len(allSequences))
    #add all sequences to table
    con = db_connect()
    cur = con.cursor()
    sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
    for seq in allSequences:
        print (seq)
        cur.execute(sequence_sql, (seq,))
    con.commit()


def getSequences():
    #return the contents of the sequence table
    conn = db_connect()
    cur = conn.cursor()
    with conn:
        cur.execute("SELECT * FROM proteins")
        return cur.fetchall()

def createMappingsHitPredict():
    #HitPredict
    conn = db_connect()
    cur = conn.cursor()
    mappings, table = extractRelevantInfoHitPredict()
    print (table.shape)
    table = table.drop_duplicates(subset = "HitPredictID")
    print (table.shape)
    match = 0
    notMatch = 0
    duplicates = 0
    for ind, row in table.iterrows():
        if ind % 100 == 0:
            print (ind)
        seq = row['Sequence']
        #retrieve myid for the matching sequence
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            #if one matching protein is found, check if the database, database_name, id pair already exists
            #print ("MATCHED")
            #id (generates automatically), database = str, database_name = str, my_id from proteins
            #print (fetched)
            my_id, _, _ = fetched[0]
            #see if already in the database
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            search = cur.execute(checkIfExists, ("HitPredict", row["HitPredictID"], my_id))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("HitPredict", row['HitPredictID'], my_id))
                match += 1
                conn.commit() #commit as it goes through, or else duplicates may be added
            else:
                print ("Already in id table")
                duplicates += 1
        else:
            #print ("OVER MATCHED")
            notMatch += 1
    print ("mathced: ", match)
    print ('duplicates: ', duplicates)
    print ("not matched: ", notMatch)
    #conn.commit()

def createMappingsHIPPIE():
    #HitPredict
    conn = db_connect()
    cur = conn.cursor()
    mappings, table = extractRelevantInfoHIPPIE()
    print (table.shape)
    #print (table.iloc[0])
    table = table.drop_duplicates(subset = "HIPPIENames")
    print (table.columns)
    #['HIPPIENames', 'Entry', 'Entry name', 'Status', 'Protein names',
    #'Cross-reference (GeneID)', 'Gene names', 'Length', 'Sequence',
    #'yourlist:M20200707216DA2B77BFBD2E6699CA9B6D1C41EB2968714N']
    print (table.shape)
    match = 0
    notMatch = 0
    duplicates = 0
    for ind, row in table.iterrows():
        if ind % 100 == 0:
            print (ind)
        seq = row['Sequence']
        #retrieve myid
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            #print ("MATCHED")
            #id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _,_ = fetched[0]
            #check if exists:
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("HIPPIE", row['HIPPIENames'], my_id))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("HIPPIE", row['HIPPIENames'], my_id))
                match += 1
                conn.commit() #commit as it goes through, or else duplicates may be added
            else:
                print ("Already in id table")
                duplicates += 1
        else:
            #print ("OVER MATCHED")
            notMatch += 1
    print ("mathced: ", match)
    print ("not matched: ", notMatch)
    print('duplicates: ', duplicates)


def createMappingsHURI():
    conn = db_connect()
    cur = conn.cursor()
    mappings, _ = openCannonicalAndIsoformFastas()
    match = 0
    notMatch = 0
    duplicates = 0
    for thing in mappings:
        seq = thing[1]
        #retrieve myid
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            #print ("MATCHED")
            #id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _, _ = fetched[0]
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("HuRI", thing[0], my_id))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("HuRI", thing[0], my_id))
                match += 1
                conn.commit()
            else:
                print ("Already in id table")
                duplicates += 1
        else:
            #print ("OVER MATCHED")
            notMatch += 1
    print ("mathced: ", match)
    print ("not matched: ", notMatch)
    print('duplicates: ', duplicates)

def addMissingSeqsHuRI():
    conn = db_connect()
    cur = conn.cursor()
    mappings, _ = openCannonicalAndIsoformFastas()
    match = 0
    notMatch = 0
    for thing in mappings:
        seq = thing[1]
        #retrieve myid
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 0:
            #print ("MATCHED")
            #id (generates automatically), database = str, database_name = str, my_id from proteins
            sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
            cur.execute(sequence_sql, (seq,))
            conn.commit()


def addEnsemblHuRIMappings():
    # add my ids
    con = db_connect()
    cur = con.cursor()
    match = 0
    notMatch = 0
    duplicates = 0
    hUniprotEnsembl = pd.read_csv("./HuRI/ensemblUniprotHuRI2.tab", sep="\t")
    eProteins = hUniprotEnsembl["Sequence"].to_list()
    # print (hUniprotEnsembl)
    hEnsembl = pd.read_csv("./HuRI/ensemblToLocation.csv")
    uniprot = hEnsembl[~(hEnsembl["Uniprot"] == "DELETED") & ~(hEnsembl["Uniprot"] == "FOLDER")]
    # print (hEnsembl)
    # merge two tables
    newDataFrame = hEnsembl.merge(hUniprotEnsembl, left_on="Uniprot",
                                  right_on="yourlist:M202007188471C63D39733769F8E060B506551E1299BFD2J")
    print(newDataFrame)
    print(newDataFrame.columns)
    # newDataFrame.to_csv("ensmeblIdsSeqs.csv")
    print(newDataFrame.shape)
    for ind, row in newDataFrame.iterrows():
        # check if each seq is added, if it is update my ids
        ensID = row["Ensembl"]
        seq = row["Sequence"]
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            # add id to the my_ids table
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _, _ = fetched[0]
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("HuRI", ensID, my_id))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                match += 1
                con.commit()
            else:
                print("Already in id table")
                duplicates += 1
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print("ensmebl uniprot matches: (out of 32)")
    print("matches: ", match)
    print("not matches: ", notMatch)
    print("duplicates: ", duplicates)
    match = 0
    notMatch = 0
    duplicates = 0
    for fasta in os.listdir("./HuRI/ensemblSequences/"):
        with open("./HuRI/ensemblSequences/" + fasta, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                pro = str(record.seq)
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                if len(fetched) == 1:
                    # add id to the my_ids table
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    my_id, _, _ = fetched[0]
                    id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                    ensID = fasta.replace(".fasta", "")
                    print(fasta, " ", ensID)
                    checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
                    cur.execute(checkIfExists, ("HuRI", ensID, my_id))
                    search = cur.fetchall()
                    if len(search) == 0:
                        id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                        cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                        match += 1
                        con.commit()
                    else:
                        print("Already in id table")
                        duplicates += 1
                else:
                    # print ("OVER MATCHED")
                    notMatch += 1
    print("folder proteins (18 total):")
    print("matches: ", match)
    print("not matches: ", notMatch)
    print("enter to commit to db ")
    input()
    con.commit()

def addHuRIEnsemblPros():
    con = db_connect()
    cur = con.cursor()
    hEnsembl = pd.read_csv("./HuRI/ensemblToLocation.csv")
    uniprot = hEnsembl[~(hEnsembl["Uniprot"] == "DELETED") & ~(hEnsembl["Uniprot"] == "FOLDER")]
    folderOnes = hEnsembl[hEnsembl["Uniprot"] == "FOLDER"]
    print(uniprot)
    print(folderOnes)
    hUniprotEnsembl = pd.read_csv("./HuRI/ensemblUniprotHuRI2.tab", sep="\t")
    eProteins = hUniprotEnsembl["Sequence"].to_list()
    # add ensembl proteins not in DB
    for pro in eProteins:
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
        fetched = cur.fetchall()
        if len(fetched) == 0:
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
            cur.execute(sequence_sql, (pro,))
            con.commit()

    # add ids c
    # folder ensembl
    for fasta in os.listdir("./HuRI/ensemblSequences/"):
        with open("./HuRI/ensemblSequences/" + fasta, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                pro = str(record.seq)
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                if len(fetched) == 0:
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                    cur.execute(sequence_sql, (pro,))
                    con.commit()
    # add my ids
    match = 0
    notMatch = 0
    duplicates= 0
    hUniprotEnsembl = pd.read_csv("./HuRI/ensemblUniprotHuRI2.tab", sep="\t")
    eProteins = hUniprotEnsembl["Sequence"].to_list()
    # print (hUniprotEnsembl)
    hEnsembl = pd.read_csv("./HuRI/ensemblToLocation.csv")
    uniprot = hEnsembl[~(hEnsembl["Uniprot"] == "DELETED") & ~(hEnsembl["Uniprot"] == "FOLDER")]
    # print (hEnsembl)
    # merge two tables
    newDataFrame = hEnsembl.merge(hUniprotEnsembl, left_on="Uniprot",
                                  right_on="yourlist:M202007188471C63D39733769F8E060B506551E1299BFD2J")
    print(newDataFrame)
    print(newDataFrame.columns)
    # newDataFrame.to_csv("ensmeblIdsSeqs.csv")
    print(newDataFrame.shape)
    for ind, row in newDataFrame.iterrows():
        # check if each seq is added, if it is update my ids
        ensID = row["Ensembl"]
        seq = row["Sequence"]
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (seq,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            # add id to the my_ids table
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _ = fetched[0]
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("HuRI", ensID, my_id))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                match += 1
                con.commit()
            else:
                print("Already in id table")
                duplicates += 1
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print ("ensmebl uniprot matches: (out of 32)")
    print("matches: ", match)
    print("not matches: ", notMatch)
    print ("duplicates: ", duplicates)
    print("enter to commit to db ")

    match = 0
    notMatch = 0
    for fasta in os.listdir("./HuRI/ensemblSequences/"):
        with open("./HuRI/ensemblSequences/" + fasta, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                pro = str(record.seq)
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                if len(fetched) == 0:
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                    cur.execute(sequence_sql, (pro,))
                elif len(fetched) == 1:
                    # add id to the my_ids table
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    my_id, _ = fetched[0]
                    id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                    ensID = fasta.replace(".fasta", "")
                    print(fasta, " ", ensID)

                    cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                    match += 1
                else:
                    # print ("OVER MATCHED")
                    notMatch += 1
    print ("folder proteins (18 total):")
    print("matches: ", match)
    print("not matches: ", notMatch)
    print ("enter to commit to db ")
    input()
    con.commit()

def addNegatomeToSqlDatbase():
    con = db_connect()
    cur = con.cursor()
    match = 0
    notMatch = 0
    negatome = pd.read_csv("./Negatome/humanSeqNegatome.tab", sep="\t")
    print(negatome.shape)
    humanNegatome = negatome[negatome["Organism ID"] == 9606]  # only human sequences
    print(humanNegatome.shape)
    addedNew = 0
    others = 0
    duplicates = 0
    for pro in humanNegatome["Sequence"].to_list():
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
        fetched = cur.fetchall()
        if len(fetched) == 0:
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
            cur.execute(sequence_sql, (pro,))
            addedNew += 1
        else:
            others += 1
    print("added: ", addedNew)
    print("already in db: ", others)
    con.commit()
    # connect negatome ids
    for ind, row in humanNegatome.iterrows():
        negId = row["NegatomeIds"]
        pro = row["Sequence"]
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            # add id to the my_ids table
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _, _ = fetched[0]

            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("Negatome", negId, my_id))
            search = cur.fetchall()
            if len(search) == 0:
                print (row)
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("Negatome", negId, my_id))
                match += 1
                con.commit()
            else:
                print("Already in id table")
                duplicates += 1
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print("matches: ", match)
    print("not matches: ", notMatch)
    print ('duplicates: ', duplicates)
    print("enter to commit to db ")
    input()
    con.commit()


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

def proteinsInDBLengthDistribution():
    #convertes the SQL table to a fasta with my_ids as the name per sequence to run 40% CD-HIT on
    distinctProQuery = "SELECT DISTINCT sequence FROM proteins"
    conn = db_connect()
    cur = conn.cursor()
    query = cur.execute(distinctProQuery)
    fetched = query.fetchall()
    #confirm that total size of db is same as unique number of seqs
    allQuery = "SELECT * FROM proteins"
    allQuery = cur.execute(allQuery)
    allFetched = allQuery.fetchall()
    lens = []
    for fetch in allFetched:
        lens.append(len(fetch[1]))
        #update table
        updateQuery = "UPDATE proteins set length=" + str(len(fetch[1])) + " where id=" + str(fetch[0])
        cur.execute(updateQuery)
    #distribution of lengths
    conn.commit()
    #sns.violinplot(x = lens)
    #plt.show()
    #distribution 2
    #sns.distplot(lens)
    #plt.show()
    #sns.distplot(lens, kde = False, rug = True)
    #plt.show()
    print ("avg")
    print (np.mean(lens))
    print ("0.5")
    print (np.percentile(lens, 5))
    print ("95")
    print (np.percentile(lens, 95))


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

def addHitPredictInteractions():
    #warning: does not check if interaction is already present,a ssumes db is empty when starting
    print ("WARNING: does not check if interaction is already present in db.  continue?")
    input()
    #fill in interaction columns for HitPredict w/ HP score
    hp = openHitPredict()
    conn = db_connect()
    cur = conn.cursor()
    match = 0
    notMatch = 0
    duplicates = 0
    updated = 0
    for ind, row in hp.iterrows():
        if ind%100 == 0:
            print ("matches: ", match,  " not matched: ", notMatch, " number duplicated pairs: ", duplicates)
        id_a = row["UniprotID_A"]
        id_b = row["UniprotID_B"]
        score = row["Confidence score"]
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HitPredict", id_a))
        fetched_a = cur.fetchall()
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HitPredict", id_b))
        fetched_b = cur.fetchall()
        if len(fetched_a) == 1 and len(fetched_b) == 1:
            #both are found
            #sort so that the pairs are ordered in the interactor db
            ids = []
            ids.append(fetched_a[0][3])
            ids.append(fetched_b[0][3])
            ids.sort()
            isHighConf = isHighConfidenceHitPredict(score)
            #see if already in table
            checkIfExists = "SELECT * FROM interaction where protein_1=? AND protein_2=? AND hitpredict_score=? AND hitpredict_high_confidence=?"
            #checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, (ids[0], ids[1], score, isHighConf))
            search = cur.fetchall()
            if len(search) == 0:
                id_map_sql = "INSERT INTO interaction (protein_1, protein_2, hitpredict_score, hitpredict_high_confidence) VALUES (?,?,?, ?)"
                cur.execute(id_map_sql, (ids[0], ids[1], score, isHighConf))
                match += 1
                conn.commit()
            else:
                duplicates += 1
        else:
            print ("FOUND: ", len(fetched_a), " protein ids interactor a ", len(fetched_b), " protein ids interactor b")
            notMatch += 1
    print ("total: ", hp.shape)
    print ("matches: ", match)
    print (" not matches: ", notMatch)


def addHIPPIEInteractions():
    #fill in interaction columns for HitPredict w/ HP score
    hp = openHIPPIE()
    print ("Total: ", hp.shape)
    conn = db_connect()
    cur = conn.cursor()
    match = 0
    updated = 0
    notMatch = 0
    duplicatePairs = []
    multIDMatches = []
    for ind, row in hp.iterrows():
        if ind%1000 == 0:
            print (ind, " out of ", hp.shape[0])
            print ("matches: ", match, " updates: ", updated, " not matched: ", notMatch, " number duplicated pairs: ", len(duplicatePairs))
        id_a = row["Uniprot Name A"]
        id_b = row["Uniprot Name B"]
        score = row["Confidence Value"]
        highConf = isHighConfidenceHIPPIE(score)
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HIPPIE", id_a))
        fetched_a = cur.fetchall()
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("HIPPIE", id_b))
        fetched_b = cur.fetchall()
        if len(fetched_a) == 1 and len(fetched_b) == 1:
            #both are found
            #sort so that the pairs are ordered in the interactor db
            ids = []
            ids.append(fetched_a[0][3])
            ids.append(fetched_b[0][3])
            ids.sort()
            #see if the combo is already in interactions
            search_interactions = "SELECT * FROM interaction WHERE protein_1=? AND protein_2=?"
            cur.execute(search_interactions, (ids[0], ids[1]))
            interaction_search = cur.fetchall()
            if len(interaction_search) == 0:
                #add new interaction
                id_map_sql = "INSERT INTO interaction (protein_1, protein_2, hippie_score, hippie_high_confidence) VALUES (?,?,?,?)"
                cur.execute(id_map_sql, (ids[0], ids[1], score, highConf))
                match += 1
                conn.commit()
            elif len(interaction_search) == 1:
                #add a HIPPIE score
                updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                #updateQuery = "UPDATE interaction set hippie_score=" + str(score) = ", hippie_high_confidence= " + str(highConf) + " where id=" + str(interaction_search[0][0])
                cur.execute(updateQuery, (score, highConf, interaction_search[0][0]))
                updated += 1
                conn.commit()
            else:
                print ("ERROR! matched more than 1 interaction......database contains duplicates....")
                duplicatePairs.append((ids[0],ids[1]))
                for matchFound in interaction_search:
                    #add HIPPIE score to each duplicate
                    updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                    #updateQuery = "UPDATE interaction set hippie_score=" + str(score) + " where id=" + str(matchFound[0])
                    #cur.execute(updateQuery)
                    cur.execute(updateQuery, (score, highConf, matchFound[0]))
                    conn.commit()
        else:
            print ("FOUND: ", len(fetched_a), " protein ids interactor a ", len(fetched_b), " protein ids interactor b")
            notMatch += 1
            if len(fetched_a) > 1:
                print ("matches mult. ids....")
                multIDMatches.append(id_a)
            elif len(fetched_b) > 1:
                print ("matches mult ids......")
                multIDMatches.append(id_b)
    print ("total: ", hp.shape)
    print ("matches: ", match)
    print (" not matches: ", notMatch)
    print ("number faulty uniprot ids....", len(multIDMatches))
    print (multIDMatches)
    return duplicatePairs


def addHuRIInteractions():
    #adds HuRI interactions to DB
    huScored = extractOnlyScoredPPIsHuRI2019()
    print (huScored.shape)
    print (huScored.iloc[0])
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
        score = row["Confidence Score"].split(":")[1] #need to split off junk "author score: #"
        if score == "-":
            print ("Error, not scored....")
        else:
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
                    id_map_sql = "INSERT INTO interaction (protein_1, protein_2, huri_score) VALUES (?,?,?)"
                    cur.execute(id_map_sql, (ids[0], ids[1], score))
                    match += 1
                    conn.commit()
                elif len(interaction_search) == 1:
                    # add a HuRI score
                    #updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                    updateQuery = "UPDATE interaction set huri_score=? WHERE id=?"
                    cur.execute(updateQuery, (score, interaction_search[0][0]))
                    updated += 1
                    conn.commit()
                else:
                    print("ERROR! matched more than 1 interaction......database contains duplicates....")
                    duplicatePairs.append((ids[0], ids[1]))
                    for matchFound in interaction_search:
                        # add HuRI score to each duplicate
                        updateQuery = "UPDATE interaction set huri_score=? WHERE id=?"
                        cur.execute(updateQuery, (score, matchFound[0]))
                        updated += 1
                        conn.commit()
                        updated += 1
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

def addNegatomeInteractions():
    #adds interactions which are marked negative in 2.0
    #if interaction already exists, sets to False
    rows = ["protein1", "protein2"]
    negatome = pd.read_csv("./Negatome/Negatome2_combinedStringent.txt", sep="\t", names = rows)
    print (negatome)
    conn = db_connect()
    cur = conn.cursor()
    match = 0
    updated = 0
    notMatch = 0
    duplicatePairs = []
    multIDMatches = []
    for ind, row in negatome.iterrows():
        if ind % 1000 == 0:
            print(ind, " out of ", negatome.shape[0])
            print("matches: ", match, " updates: ", updated, " not matched: ", notMatch, " number duplicated pairs: ",
                  len(duplicatePairs), " multi  matched HuRI IDs: ", len(multIDMatches))
        id_a = row["protein1"]
        id_b = row["protein2"]
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("Negatome", id_a))
        fetched_a = cur.fetchall()
        cur.execute("SELECT * FROM id_map WHERE database=? AND database_name=?", ("Negatome", id_b))
        fetched_b = cur.fetchall()
        if len(fetched_a) == 1 and len(fetched_b) == 1:
            # both are found
            # sort so that the pairs are ordered in the interactor db
            ids = []
            ids.append(fetched_a[0][3])
            ids.append(fetched_b[0][3])
            ids.sort()
            # print (ids)
            # print (score)
            # see if the combo is already in interactions
            search_interactions = "SELECT * FROM interaction WHERE protein_1=? AND protein_2=?"
            cur.execute(search_interactions, (ids[0], ids[1]))
            interaction_search = cur.fetchall()
            if len(interaction_search) == 0:
                # add new interaction
                id_map_sql = "INSERT INTO interaction (protein_1, protein_2, negatome_interaction) VALUES (?,?,?)"
                cur.execute(id_map_sql, (ids[0], ids[1], True))
                match += 1
                conn.commit()
            elif len(interaction_search) == 1:
                # add a HuRI score
                # updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                updateQuery = "UPDATE interaction set negatome_interaction=? WHERE id=?"
                cur.execute(updateQuery, (True, interaction_search[0][0]))
                updated += 1
                conn.commit()
            else:
                print("ERROR! matched more than 1 interaction......database contains duplicates....")
                duplicatePairs.append((ids[0], ids[1]))
                for matchFound in interaction_search:
                    # add HuRI score to each duplicate
                    updateQuery = "UPDATE interaction set negatome_interaction=? WHERE id=?"
                    cur.execute(updateQuery, (True, matchFound[0]))
                    updated += 1
                    conn.commit()
                    updated += 1
        else:
            notMatch += 1
            print ("FOUND: ", len(fetched_a), " ", len(fetched_b))
    print("matches: ", match)
    print(" not matches: ", notMatch)
    print ("updated: ", updated)

def addNegatomeIDsFixed():
    negatome = pd.read_csv("./Negatome/humanSeqNegatome.tab", sep="\t")
    print(negatome.shape)
    humanNegatome = negatome[negatome["Organism ID"] == 9606]  # only human sequences
    print(humanNegatome)
    expanded = {'negatomeID': [], 'sequence': []}
    for ind, row in humanNegatome.iterrows():
        ids = row['NegatomeIds']
        seq = row['Sequence']
        trySplit = ids.split(",")
        if len(trySplit) == 1:
            expanded['negatomeID'].append(trySplit[0])
            expanded['sequence'].append(seq)
        else:
            for i in trySplit:
                expanded['negatomeID'].append(i)
                expanded['sequence'].append(seq)
    asDF = pd.DataFrame(expanded)
    asDF = asDF.drop_duplicates()
    badNames = ['P0CG05', 'P62988']
    asDF = asDF[~(asDF['negatomeID'].isin(badNames))]
    print(asDF.negatomeID.value_counts())
    print(asDF)
    #make sure the sequences are all added
    con = db_connect()
    cur = con.cursor()
    match= 0
    notMatch = 0
    duplicates = 0
    for ind, row in asDF.iterrows():
        negId = row['negatomeID']
        pro = row['sequence']
        cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
        fetched = cur.fetchall()
        if len(fetched) == 1:
            # add id to the my_ids table
            # print ("MATCHED")
            # id (generates automatically), database = str, database_name = str, my_id from proteins
            my_id, _, _ = fetched[0]
            checkIfExists = "SELECT * FROM id_map where database =? AND database_name=? AND my_id=?"
            cur.execute(checkIfExists, ("Negatome", negId, my_id))
            search = cur.fetchall()
            if len(search) == 0:
                print(row)
                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                cur.execute(id_map_sql, ("Negatome", negId, my_id))
                match += 1
                con.commit()
            else:
                print("Already in id table")
                duplicates += 1
        else:
            # print ("OVER MATCHED")
            notMatch += 1
    print("matches: ", match)
    print("not matches: ", notMatch)
    print('duplicates: ', duplicates)
    print("enter to commit to db ")


def fixCommas():
    #if commas are present in the id map table, removes them and makes new entries
    query = "SELECT * FROM id_map WHERE database_name LIKE '%,%'"
    con = db_connect()
    cur = con.cursor()
    cur.execute(query)
    results = cur.fetchall()
    print(results)
    badIds = []
    for result in results:
        badIds.append(result[0])
        db = result[1]
        toSplit = result[2]
        proID = result[3]
        splits = toSplit.split(",")
        for thing in splits:
            query = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
            cur.execute(query, (db, thing, proID))
            con.commit()
    for badID in badIds:
        query = "DELETE FROM id_map WHERE id=?"
        cur.execute(query, (badID,))
        con.commit()

def getAllProteinIDsInDBHuRI():
    #return set of all protein ids in the DB so far
    conn = db_connect()
    cur = conn.cursor()
    query = "SELECT * FROM id_map"
    cur.execute(query)
    retreived = cur.fetchall()
    return set([x[2] for x in retreived if x[1] == "HuRI"])

def getAllProteinIDsInDB():
    #return set of all protein ids in the DB so far
    conn = db_connect()
    cur = conn.cursor()
    query = "SELECT * FROM id_map"
    cur.execute(query)
    retreived = cur.fetchall()
    return set([x[2] for x in retreived])

#adding more HuRI interactions (including negatives)
def openTotalHuRIProteins():
    #get HuRI proteins not in DB
    openHuRI = openMitab27File("./HuRI/HI-union.psi")
    #print (openHuRI)
    idsA = [x.split(':')[1] for x in openHuRI['Unique ID A'].to_list()]
    idsA = list(set(idsA))
    idsB = [x.split(':')[1] for x in openHuRI['Unique ID B'].to_list()]
    idsB = list(set(idsB))
    allIds = set(idsA + idsB)
    print (len(allIds))
    #figure out what's missing
    allDBPros = getAllProteinIDsInDB()
    notInDB2 = allIds.difference(allDBPros)
    inDB = allIds.intersection(allDBPros)
    print ("those in DB: ", len(inDB))
    print("not found in DB: ", len(notInDB2))
    return notInDB2, inDB

def updateHuRIMappings():
    #with those in DB, check to make sure the mappings are all there
    openHuRI = openMitab27File("./HuRI/HI-union.psi")
    # print (openHuRI)
    idsA = [x.split(':')[1] for x in openHuRI['Unique ID A'].to_list()]
    idsA = list(set(idsA))
    idsB = [x.split(':')[1] for x in openHuRI['Unique ID B'].to_list()]
    idsB = list(set(idsB))
    allIDs = set(idsA + idsB)
    thoseInHuRI = getAllProteinIDsInDBHuRI()
    HuRINotPresent = allIDs.difference(thoseInHuRI)
    print (len(HuRINotPresent)) #1441 not in DB with HuRI mark
    allInDB = getAllProteinIDsInDB()
    thoseInDBNotMarkedHuRI = HuRINotPresent.intersection(allInDB)
    print (len(thoseInDBNotMarkedHuRI)) #396 of these, add to id_map
    conn = db_connect()
    cur = conn.cursor()
    for id in list(thoseInDBNotMarkedHuRI):
        #print (id)
        #get id associated with that protein
        query = "SELECT * FROM id_map where database_name=?"
        cur.execute(query, (id,))
        fetched = cur.fetchall()
        #print (fetched)
        proID = fetched[0][3]
        #print ("pro id: ", proID)
        #check exists
        query = "SELECT * FROM id_map where database = ? AND database_name=?"
        cur.execute(query, ("HuRI", id))
        fetched = cur.fetchall()
        #print (fetched)
        if fetched:
            print ("ERROR!")
        else:
            #update DB
            query = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
            cur.execute(query, ("HuRI", id, proID))
    conn.commit()




def openCannonicalAndIsoformFastas(cannonical, noncannonical):
    #takes a list of isoform ids and pulls the sequence from teh Uniprot alternative isoforms fasta file
    #turns into dictionary with key:sequence entries
    #cannonical, noncannonical = getListAllAlternativeIsoforms()
    keeperSeqs = []
    foundIds = [] #mappning of tuples from HuRI to uniprot ids
    idListIntermed  = []
    #print (len(noncannonical))
    with open("./HuRI/cannonicalAndIsoformsHuRI.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            if actualID in noncannonical:
                keeperSeqs.append(str(record.seq))
                foundIds.append((actualID, str(record.seq)))
                idListIntermed.append(actualID)
            elif actualID in cannonical:
                keeperSeqs.append(str(record.seq))
                foundIds.append((actualID, str(record.seq)))
                idListIntermed.append(actualID)
    #print ("found: ", len(keeperSeqs))
    #figure out what is not there
    foundSet = set(idListIntermed)
    notFound = list(set(noncannonical).difference(foundSet))
    notFoundCannon = list(set(cannonical).difference(foundSet))
    print ("len not found cannonical: ", len(notFoundCannon))
    print("len not found noncannonical: ", len(notFound))
    justOnesNotFoundCannon = [x for x in notFoundCannon if "-1" in x]
    print ("len not found cannonical w/ -1: ", len(justOnesNotFoundCannon))
    #print (len(notFound))
    with open("./HuRI/isoformFastas.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            #print ("id: ", actualID)
            for missing in notFound:
                stub = missing.split("-")[0]
                if stub == actualID:
                    keeperSeqs.append(str(record.seq))
                    foundIds.append((missing, str(record.seq)))
                    idListIntermed.append(missing)
    #print("found: ", len(keeperSeqs))
    #print ("found: ", len(set(keeperSeqs)))
    #print ("ids: ", len(foundIds))
    #getting -1 isoforms not found
    with open("./HuRI/cannonicalAndIsoformsHuRI.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            for missing in justOnesNotFoundCannon:
                stub = missing.split("-")[0]
                if stub == actualID:
                    #looking for -1 isoforms left over which are not
                    keeperSeqs.append(str(record.seq))
                    foundIds.append((missing, str(record.seq)))
                    idListIntermed.append(missing)
    foundSet = set(foundIds)
    notFound = list(set(noncannonical).difference(foundSet))
    print ("total ids: ", len(cannonical) + len(noncannonical))
    print ("Found: ", len(foundIds))
    print ("found seqs: ", len(keeperSeqs))
    totalIds = set(cannonical + noncannonical)
    notFound = list(totalIds.difference(set(idListIntermed)))
    print ("total not found: ", len(notFound))
    print("total ids: ", len(cannonical) + len(noncannonical))
    print("Found: ", len(set(idListIntermed)))
    print("found seqs: ", len(set(keeperSeqs)))
    return foundIds, keeperSeqs

def nametype(x):
    #returns if HuRI or Ensembl name
    if "." in x:
        return "Ensembl"
    else:
        return "HuRI"

def isoformID(x):
    if "." not in x:
        if "-" in x:
            return int(x.split("-")[1])
        else:
            return 1
    else:
        print ("no isofrm ")
        return None


def addHuRIInteractionsV2():
    #adds HuRI interactions to DB
    huScored = openMitab27File("./HuRI/HI-union.psi")
    print (huScored.shape)
    print (huScored.iloc[0])
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
        score = row["Confidence Score"].split(":")[1] #need to split off junk "author score: #"
        if score == "-":
            score = -1
            #adding all HuRI
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
                    id_map_sql = "INSERT INTO interaction (protein_1, protein_2, huri_score) VALUES (?,?,?)"
                    cur.execute(id_map_sql, (ids[0], ids[1], -1))
                    match += 1
                    conn.commit()
                elif len(interaction_search) == 1:
                    # add a HuRI score
                    #updateQuery = "UPDATE interaction SET hippie_score=?, hippie_high_confidence=? WHERE id=?"
                    updateQuery = "UPDATE interaction set huri_score=? WHERE id=?"
                    cur.execute(updateQuery, (score, interaction_search[0][0]))
                    updated += 1
                    conn.commit()
                else:
                    print("ERROR! matched more than 1 interaction......database contains duplicates....")
                    duplicatePairs.append((ids[0], ids[1]))
                    for matchFound in interaction_search:
                        # add HuRI score to each duplicate
                        updateQuery = "UPDATE interaction set huri_score=? WHERE id=?"
                        cur.execute(updateQuery, (score, matchFound[0]))
                        updated += 1
                        conn.commit()
                        updated += 1
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

def updateHuRIEnsemblSeqs():
    con = db_connect()
    cur = con.cursor()
    match = 0
    notMatch = 0
    prosNew = 0
    prosNotNew  = 0
    alreadyFound= 0
    for fasta in os.listdir("./HuRI/ensemblSequences/"):
        with open("./HuRI/ensemblSequences/" + fasta, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                pro = str(record.seq)
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                print ("found protein: ", fetched)
                if len(fetched) == 0:
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                    cur.execute(sequence_sql, (pro,))
                    print ("entering new protein ")
                    prosNew += 1
                elif len(fetched) == 1:
                    prosNotNew += 1
                    # add id to the my_ids table
                    # print ("MATCHED")
                    # id (generates automatically), database = str, database_name = str, my_id from proteins
                    my_id, _, _ = fetched[0]
                    # check if it already exists
                    checkQuery = "SELECT * FROM id_map WHERE database=? and database_name=? and my_id=?"
                    ensID = fasta.replace(".fasta", "")
                    print(fasta, " ", ensID)
                    cur.execute(checkQuery, ("HuRI", ensID, my_id))
                    if len(cur.fetchall()) == 0:
                        id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                        cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                        match += 1
                    else:
                        print ("already found")
                        alreadyFound += 1
                else:
                    # print ("OVER MATCHED")
                    notMatch += 1
    print("folder proteins (18 total):")
    print("matches: ", match)
    print("not matches: ", notMatch)
    print ("new pros: ", prosNew)
    print ("found pros: ", prosNotNew)
    print("enter to commit to db ")
    print (" already in id_map: ", alreadyFound)
    input()
    con.commit()

def addingHuRIFastaPart2():
    #adding uniprot seqs in text file which were not present in the DB already
    #in fasta in folder
    #all names to try and get from the fasta (some failed again)
    toAdd = filesToGrab = open("remainingHuRIUniprot.txt", "r")
    toAdd = [x.strip() for x in toAdd.readlines()]
    con = db_connect()
    cur = con.cursor()
    prosNew = 0
    prosNotNew = 0
    allTheFastaNames = []
    seqDict = {'name':[], 'seq':[]}
    with open("./HuRI/uniprotHuRIRound2.fasta", mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            actualID = record.id.split("|")[1]
            allTheFastaNames.append(actualID)
            seqDict['name'].append(actualID)
            seqDict['seq'].append(str(record.seq))
            '''
            if actualID in toAdd:
                #check again not in DB
                pro = str(record.seq)
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                #print("found protein: ", fetched)
                if len(fetched) == 0:
                        #add seq to DB
                        sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                        cur.execute(sequence_sql, (pro,))
                        #print("entering new protein ")
                        prosNew += 1
                else:
                    #WTF???
                    prosNotNew += 1
            '''
    toAddWithOne = [x for x in toAdd if "-1" in x]
    otherToAdd = [x for x in toAdd if "-1" not in x]
    print ("total in to add: ", len(toAdd), " with -1: ", len(toAddWithOne), " other in to add ", len(otherToAdd))
    allEndInOneFasta = [x for x in allTheFastaNames if "-1" in x]
    otherInFasta = [x for x in allTheFastaNames if "-1" not in x]
    print("total in to add: ", len(allTheFastaNames), " with -1: ", len(allEndInOneFasta), " other in to add ", len(otherInFasta))

    #figure out which proteins lacking isoform 1 and add it manually
    fastaBaseNames = [x.split("-")[1] for x in allTheFastaNames if "-" in x]
    otherBaseNames = [x for x in allTheFastaNames if "-" not in x]
    totalBaseNames = set(list(set(fastaBaseNames)) + list(set(otherBaseNames)))
    print ('unique bases ', len(totalBaseNames))
    addedOneReps = []
    for x in totalBaseNames:
        if x + "-1" not in allTheFastaNames:
            #track this for x-1 option
            addedOneReps.append(x)
    print (len(addedOneReps))
    #combine toAdd and addedOneReps
    #comboReps = list(set(toAdd + addedOneReps))
    seqDF = pd.DataFrame(seqDict)
    seqsFound = 0
    seqsNotFound = 0
    match = 0
    alreadyFound = 0
    for uniprotName in toAdd:
        print (uniprotName)
        if "-1" in uniprotName:
            if uniprotName.split("-")[0] in addedOneReps:
                matching = seqDF[seqDF.name == uniprotName.split("-")[0]]
                if matching.shape[0] != 0:
                    pro = matching['seq'].to_list()[0]
                    #print("FOUND SEQ")
                    #print(matching)
                    seqsFound += 1
                    #check if in DB
                    cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                    fetched = cur.fetchall()
                    # print("found protein: ", fetched)
                    if len(fetched) == 0:
                        # add seq to DB
                        sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                        cur.execute(sequence_sql, (pro,))
                        # print("entering new protein ")
                        prosNew += 1
                    else:
                        # WTF???
                        prosNotNew += 1
                        my_id, _, _ = fetched[0]
                        # check if it already exists
                        checkQuery = "SELECT * FROM id_map WHERE database=? and database_name=? and my_id=?"
                        ensID = uniprotName
                        cur.execute(checkQuery, ("HuRI", ensID, my_id))
                        if len(cur.fetchall()) == 0:
                            id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                            cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                            match += 1
                        else:
                            print("already found")
                            alreadyFound += 1
                else:
                    #print("NOT FOUND")
                    seqsNotFound += 1
            else:
                matching = seqDF[seqDF.name == uniprotName]
                if matching.shape[0] != 0:
                    #print("FOUND SEQ")
                    #print(matching)
                    seqsFound += 1
                    pro = matching['seq'].to_list()[0]
                    #print("FOUND SEQ")
                    #print(matching)
                    # check if in DB
                    cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                    fetched = cur.fetchall()
                    # print("found protein: ", fetched)
                    if len(fetched) == 0:
                        # add seq to DB
                        sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                        cur.execute(sequence_sql, (pro,))
                        #print("entering new protein ")
                        prosNew += 1
                    else:
                        # WTF???
                        prosNotNew += 1
                        # check if in DB
                        cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                        fetched = cur.fetchall()
                        # print("found protein: ", fetched)
                        if len(fetched) == 0:
                            # add seq to DB
                            sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                            cur.execute(sequence_sql, (pro,))
                            # print("entering new protein ")
                            prosNew += 1
                        else:
                            # WTF???
                            prosNotNew += 1
                            my_id, _, _ = fetched[0]
                            # check if it already exists
                            checkQuery = "SELECT * FROM id_map WHERE database=? and database_name=? and my_id=?"
                            ensID = uniprotName
                            cur.execute(checkQuery, ("HuRI", ensID, my_id))
                            if len(cur.fetchall()) == 0:
                                id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                                cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                                match += 1
                            else:
                                print("already found")
                                alreadyFound += 1
                else:
                        #print("NOT FOUND")
                        seqsNotFound += 1
        else:
            #not a 1 isoform, just check normal
            matching = seqDF[seqDF.name == uniprotName]
            if matching.shape[0] != 0:
                #print("FOUND SEQ")
                #print(matching)
                seqsFound += 1
                pro = matching['seq'].to_list()[0]
                #print("FOUND SEQ")
                #print(matching)
                # check if in DB
                cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                fetched = cur.fetchall()
                # print("found protein: ", fetched)
                if len(fetched) == 0:
                    # add seq to DB
                    sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                    cur.execute(sequence_sql, (pro,))
                    # print("entering new protein ")
                    prosNew += 1
                else:
                    prosNotNew += 1
                    # check if in DB
                    cur.execute("SELECT * FROM proteins WHERE sequence=?", (pro,))
                    fetched = cur.fetchall()
                    # print("found protein: ", fetched)
                    if len(fetched) == 0:
                        # add seq to DB
                        sequence_sql = "INSERT INTO proteins (sequence) VALUES (?)"
                        cur.execute(sequence_sql, (pro,))
                        # print("entering new protein ")
                        prosNew += 1
                    else:
                        # WTF???
                        prosNotNew += 1
                        my_id, _, _ = fetched[0]
                        # check if it already exists
                        checkQuery = "SELECT * FROM id_map WHERE database=? and database_name=? and my_id=?"
                        ensID = uniprotName
                        cur.execute(checkQuery, ("HuRI", ensID, my_id))
                        if len(cur.fetchall()) == 0:
                            id_map_sql = "INSERT INTO id_map (database, database_name, my_id) VALUES (?,?,?)"
                            cur.execute(id_map_sql, ("HuRI", ensID, my_id))
                            match += 1
                        else:
                            print("already found")
                            alreadyFound += 1
            else:
                #print("NOT FOUND")
                seqsNotFound += 1
    print ("Found: ", seqsFound)
    print ("Seqs not found: ", seqsNotFound)
    print ("pros new: ", prosNew)
    print ("pros not new: ", prosNotNew)
    print ("matched: ", match)
    print (" not matched: ", alreadyFound)
    input()
    con.commit()



def addHuRIInteractionsPart2():
    #adds HuRI interactions to DB that have no score
    huScored = openMitab27File("./HuRI/HI-union.psi")
    print(huScored.columns)
    print("HURI TOTAL: ", huScored.shape[0])
    huScored = huScored[huScored['Confidence Score'] == "-"] #only those not scored already, as those are already in the DB
    print (huScored.shape)
    print (huScored.iloc[0])
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
                    updated += 1
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


#add huri negatives (?)


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    #proteinsInDBLengthDistribution()
    #update HuRI positives set
    #addHuRIInteractionsV2()
    '''
    notinDB, inDB = openTotalHuRIProteins()
    print (len(notinDB))
    for thing in notinDB:
        print (thing)
    uniprotNotinDB = [x for x in notinDB if "." not in x]
    ensemblID = [x for x in notinDB if "." in x]
    print (len(uniprotNotinDB))
    filesToGrab = open("remainingHuRIUniprot.txt", "w")
    filesToGrab.write("\n".join(uniprotNotinDB))
    filesToGrab.close()
    filesToGrab2  = open("remainingHuRIEnsembl.txt", "w")
    filesToGrab2.write("\n".join(ensemblID))
    filesToGrab2.close()
    print (len(ensemblID))
    #print ()
    '''
    #updateHuRIEnsemblSeqs()
    #
    #addHuRIInteractionsPart2()
    '''
    con = db_connect()
    cur = con.cursor()
    addNewCol = "ALTER TABLE interaction ADD COLUMN huri_negative boolean"
    cur.execute(addNewCol)
    con.commit()
    '''

    #add interactions of HuRI
    '''
    huScored = openHuRI = openMitab27File("./HuRI/HI-union.psi")
    print (huScored.columns)
    print ("HURI TOTAL: ", huScored.shape[0])
    huNoTScored = huScored[huScored['Confidence Score'] == "-"]
    print ("HuRI not scored: ", huNoTScored.shape[0])
    print ("HuRI scored: ", huScored.shape[0] - huNoTScored.shape[0])
    print (huScored['Confidence Score'].value_counts())
    alreadyDone = extractOnlyScoredPPIsHuRI2019()
    print (alreadyDone.shape[0])
    huScored = openHuRI = openMitab27File("./HuRI/HuRI.psi")
    print (huScored.shape[0])
    addHuRIInteractionsPart2()
    '''
