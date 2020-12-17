#getting Fragoza sequences and mutant sequences from NCBI

# AML

# 6.10.20

# NM to NP, pull NP as fastas
import pandas as pd
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

# set up Entrez usage
Entrez.email = "lafleur1@cs.washington.edu"
Entrez.api_key = None  # note that the 3 requests per second cut off is enforced automatically by BioPython

# keys for ncbi nucleotide
'''
handle = Entrez.einfo(db="nucleotide")
record = Entrez.read(handle)
record["DbInfo"]["Description"]

for field in record["DbInfo"]["FieldList"]:
	print("%(Name)s, %(FullName)s, %(Description)s" % field)

ALL, All Fields, All terms from all searchable fields
UID, UID, Unique number assigned to each sequence
FILT, Filter, Limits the records
WORD, Text Word, Free text associated with record
TITL, Title, Words in definition line
KYWD, Keyword, Nonstandardized terms provided by submitter
AUTH, Author, Author(s) of publication
JOUR, Journal, Journal abbreviation of publication
VOL, Volume, Volume number of publication
ISS, Issue, Issue number of publication
PAGE, Page Number, Page number(s) of publication
ORGN, Organism, Scientific and common names of organism, and all higher levels of taxonomy
ACCN, Accession, Accession number of sequence
PACC, Primary Accession, Does not include retired secondary accessions
GENE, Gene Name, Name of gene associated with sequence
PROT, Protein Name, Name of protein associated with sequence
ECNO, EC/RN Number, EC number for enzyme or CAS registry number
PDAT, Publication Date, Date sequence added to GenBank
MDAT, Modification Date, Date of last update
SUBS, Substance Name, CAS chemical name or MEDLINE Substance Name
PROP, Properties, Classification by source qualifiers and molecule type
SQID, SeqID String, String identifier for sequence
GPRJ, BioProject, BioProject
SLEN, Sequence Length, Length of sequence
FKEY, Feature key, Feature annotated on sequence
PORG, Primary Organism, Scientific and common names of primary organism, and all higher levels of taxonomy
COMP, Component Accession, Component accessions for an assembly
ASSM, Assembly, Assembly
DIV, Division, Division
STRN, Strain, Strain
ISOL, Isolate, Isolate
CULT, Cultivar, Cultivar
BRD, Breed, Breed
BIOS, BioSample, BioSample


'''


def fragozaTranscriptIDs():
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

    print(origDataset.HGVS_cDNA)
    justNucleotideUIDs = origDataset.apply(lambda row: row['HGVS_cDNA'].split(':')[0], axis=1)
    print(justNucleotideUIDs)
    # collapse to unique UIDs
    justNucleotideUIDs = justNucleotideUIDs.drop_duplicates()
    # checking no nulls exist
    print("number null rows: ", justNucleotideUIDs.isnull().sum())
    # for each NT ID, pull the nucleotide entry
    return justNucleotideUIDs


def createTranscriptToProteinKey():
    ntIDs = fragozaTranscriptIDs()
    listIDs = ntIDs.to_list()

    storageDF = {'transcript': [], 'number products': [], 'product': []}
    products = []
    numberProductsAssociatedWithTranscript = []
    # ntIDs = listIDs[0:5]
    # for each, pull and get the NP associated as a fasta
    for i, ntID in enumerate(listIDs):
        # try with one and see output
        fetched = Entrez.efetch(db='nucleotide', id=ntID, retmode='xml')
        parsed = Entrez.read(fetched)
        print("on ", i, " out of ", len(listIDs), " ", ntID)
        foundCDS = False
        for key in parsed[0].keys():
            # print (key)
            # print (parsed[0][key])
            # print ( " ")
            if key == 'GBSeq_feature-table':
                for entry in parsed[0][key]:
                    # print (entry)
                    # print (entry['GBFeature_key'])
                    if entry['GBFeature_key'] == 'CDS':
                        intervals = [x for x in entry['GBFeature_quals'] if x['GBQualifier_name'] == 'protein_id']
                        print("Found: ", len(intervals), " protiein products")
                        print(intervals)
                        numberProductsAssociatedWithTranscript.append(len(intervals))
                        justIds = [x['GBQualifier_value'] for x in intervals]
                        asString = ",".join(justIds)
                        products.append(asString)
                        foundCDS = True
                    # for thing in entry:
                    # print (thing)
                    # print (entry[thing])
        # if never found a CDS, add dummy values
        if not foundCDS:
            numberProductsAssociatedWithTranscript.append(0)
            products.append("")

    print(ntIDs)
    print(numberProductsAssociatedWithTranscript)
    print(products)
    storageDF = {'transcript': listIDs, 'number products': numberProductsAssociatedWithTranscript, 'product': products}
    storageDF = pd.DataFrame(storageDF)
    storageDF.to_csv("transcriptProteinTableFragoza.csv", index=False)


'''
GBSeq_locus
GBSeq_length
GBSeq_strandedness
GBSeq_moltype
GBSeq_topology
GBSeq_division
GBSeq_update-date
GBSeq_create-date
GBSeq_definition
GBSeq_primary-accession
GBSeq_accession-version
GBSeq_other-seqids
GBSeq_keywords
GBSeq_source
GBSeq_organism
GBSeq_taxonomy
GBSeq_references
GBSeq_comment
GBSeq_primary
GBSeq_feature-table
GBSeq_sequence
'''


def pullAllFragozaFastas():
    openedKey = pd.read_csv("transcriptProteinTableFragoza.csv")
    fastaDir = "./FragozaFastas/"
    for name in openedKey['product']:
        saveName = fastaDir + name + '.fa'
        # download FASTA and save
        fetched = Entrez.efetch(db='protein', id=name, rettype="fasta", retmode="text", )
        data = fetched.read()
        print(data)
        fetched.close()
        outHandle = open(saveName, 'w')
        outHandle.write(data)
        outHandle.close()
# save fastas to a file


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

targetIDs = origDataset['Target Entrez GeneID'].drop_duplicates()
targetSet = set(targetIDs.to_list())
interactorIDs = origDataset['Interactor Entrez GeneID'].drop_duplicates()
interactorSet = set(interactorIDs.to_list())
print(targetIDs.shape)
print(interactorIDs.shape)
print("only in target: ", len(targetSet.difference(interactorSet)))
print("only in interactor: ", len(interactorSet.difference(targetSet)))
merged = targetIDs.to_list() + interactorIDs.to_list()
print(len(set(merged)))
print("in both: ",
      len(set(merged)) - len(targetSet.difference(interactorSet)) - len(interactorSet.difference(targetSet)))

interactorIDs = open("interactorGeneIDs.txt", "w")
# printing just the interactor set to file to pull NP ids for the cannonical isoforms
for i in interactorSet.difference(targetSet):
    interactorIDs.write(i)
    interactorIDs.write("\n")
interactorIDs.close()

mappings = pd.read_csv("mapping.tab", sep='\t')
print(mappings)
counts = mappings['yourlist'].value_counts()
print(counts[counts > 1])

'''
122183    5
441521    3
1409      2
51207     2
9465      2

'''

# ddg files
# 1A22_A	1A22_B	1A22_A_CA171A	1A22_B	-12.994037 -11.984037 1.01
# wt_a wt_b mt_a mt_b score_wt score_mt diff (per line)

# seq format
# name (same as ddg) tab sequence

# in fastas: 400 + 387
# in table: all rest but above

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

print(origDataset.HGVS_cDNA)
origDataset['transcript'] = origDataset.apply(lambda row: row['HGVS_cDNA'].split(':')[0], axis=1)
# origDataset['targetBase'] = origDataset.apply(lambda row: row['UniProt'].split(',')[0], axis = 1)
origDataset['mutationAA'] = origDataset.apply(lambda row: row['UniProt'].split(',')[1], axis=1)

print(origDataset)
openedKey = pd.read_csv("transcriptProteinTableFragoza.csv")
merged = pd.merge(origDataset, openedKey, how='left', on='transcript')
print(merged)

print(merged['product'].isnull().sum())

merged['saveIDTarget'] = merged['Target Entrez GeneID'] + "_" + merged['mutationAA']
print(merged)


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


tooLong = []
seqTable = {"ID": [], "Seq": [], 'Len': []}
# saveIDTarget (mutation partner) mut id
# "Interactor Entrez GeneID" other
merged.to_csv("mergedIntermed")
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
        print("Bad mutation")
        badMutants.append(saveID)
    if len(wtSeq) > 1:
        seqTable["ID"].append(origID)
        seqTable["Seq"].append(wtSeq)
        seqTable['Len'].append(len(wtSeq))
print(merged[merged.saveIDTarget.isin(badMutants)])
print(badMutants)

# for each row go through and add interactor WT to seq table

badNamesInteractors = [122183, 441521, 1409, 51207, 9465]

# saving those interactor wt to a dataframe
uniprotInteractors = pd.read_csv("withIDsInteractorUniprot.tab", sep='\t')
print(uniprotInteractors)
# now for each not in the bad set, get the seqT
for i, row in uniprotInteractors.iterrows():
    idInt = row['yourlist']
    # print (type(idInt))
    sequence = row['Sequence']
    if idInt not in badNamesInteractors:
        seqTable['ID'].append(idInt)
        seqTable['Seq'].append(sequence)
        seqTable['Len'].append(len(sequence))

# add sequences from interactor
# manually gahtering fastas of bad ones
'''
Not mapped:
not mapped
390535 (psuedogene, nontranslated)
10597 (done)
285733 (done)
541471  (HuRI 14, not found....) 
145946	(done)
729862 (psuedogene, nontranslated)
440184  (done)
202459   (psuedogene, not found)
155060 (psuedogene, not found)
440321 (done)
401508 (psuedogene, not found)
9753 (done) 
122183    5 (done) 
441521    3 (done)
1409      2 (done)
51207     2 (done)
9465      2 (done)

'''
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
print(asDF.shape)
asDF = asDF[asDF['Len'] <= 600]
print(asDF.shape)
passedIDS = asDF['ID'].to_list()
asDF = asDF.drop(['Len'], axis=1)
asDF.to_csv("FragozaSequencesML600.tab", sep="\t", header=False, index=False)


def flipDisruption(disr):
    if disr == 1:
        return 0
    else:
        return 1


# toPredict
ppiPairs = {'wtOne': [], 'wtTwo': [], 'mtOne': [], 'mtTwo': [], 'wtDG': [], 'mtDG': [], 'delta': []}
PIPRPairs = {'v1': [], 'v2': [], 'label': [], 'type': []}  # 1 if interacts, 0 else
print("TOO LONG COUNT: ", len(tooLong))
for i, row in merged.iterrows():
    # open orig FASTA
    origID = row['Target Entrez GeneID']  # wt 1
    saveID = row['saveIDTarget']  # mt 1
    wtID2 = row["Interactor Entrez GeneID"]
    # check if not one I could find
    # print (type(wtID2))
    # passOrig = origID not in tooLong and wtID2 not in tooLong and wtID2 not in badInteractors
    # passMt = saveID not in badMutants and saveID not in tooLong and wtID2 in tooLong
    # if wtID2 not in badInteractors and saveID not in badMutants and wtID2 not in tooLong and origID not in tooLong:
    if wtID2 in passedIDS and origID in passedIDS and saveID in passedIDS:
        ppiPairs['wtOne'].append(origID)
        ppiPairs['wtTwo'].append(wtID2)
        ppiPairs['mtOne'].append(saveID)
        ppiPairs['mtTwo'].append(wtID2)
        ppiPairs['wtDG'].append(0)
        ppiPairs['mtDG'].append(0)
        disruption = row['Disruption']  # 0 if nondisrupting, 1 if disrupting
        ppiPairs['delta'].append(disruption)

        # PIPR input (max len 600)
        PIPRPairs['v1'].append(origID)
        PIPRPairs['v2'].append(wtID2)
        PIPRPairs['label'].append(1)
        PIPRPairs['type'].append("WT")
        PIPRPairs['v1'].append(saveID)
        PIPRPairs['v2'].append(wtID2)
        print("DIS WAS: ", disruption)
        print("DIS NOW: ", flipDisruption(disruption))
        PIPRPairs['label'].append(flipDisruption(disruption))
        PIPRPairs['type'].append("MT")




asDF = pd.DataFrame(ppiPairs)
print(asDF)
print(asDF.delta.value_counts())
asDF.to_csv("ddgDummyFragozaML600.tab", sep="\t", header=False, index=False)


"""
# 868 non disruptive mutants, 195 disrutpive

piprDF = pd.DataFrame(PIPRPairs)
print(piprDF)
piprDF = piprDF.drop_duplicates()
print(piprDF)

# save only WT
piprDFWT = piprDF[piprDF['type'] == 'WT']
print(piprDFWT)
piprDFWT = piprDFWT.drop(['type'], axis=1)
piprDFWT.to_csv("wtFragozaPIPR.tab", sep="\t", header=True, index=False)
# save only MT
# save only WT
piprDFWT = piprDF[piprDF['type'] == 'MT']
print(piprDFWT)
piprDFWT = piprDFWT.drop(['type'], axis=1)
piprDFWT.to_csv("mtFragozaPIPR.tab", sep="\t", header=True, index=False)
# save all
# save only WT

piprDF = piprDF.drop(['type'], axis=1)
piprDF.to_csv("bothFragozaPIPR.tab", sep="\t", header=True, index=False)
"""

#looking at how many of these sequencs alreayd in the PPI DB

#looking at intersection of fragoza sequences and ppiDB sequences:



