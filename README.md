# embeddingsPPI: Binary protein-protien interaction prediction with unsupervised learning embeddings

# Example Project

Protein-protein interactions (PPIs) are the main mechanism by which proteins carry out their functions.  A better ability predict if two proteins will interact can be 

VIDEO GOES HERE

## Introduction

Overview of protein sequences.  PPIs.  Why predicting is important. 
Mutatnt ppi stuff.

## Related Work

PIPR.  MuPIPR.  Some embeddings used.  Best performacne as of yet.

Introduce Fragoza dataset. Show poor perofrmance on Fragoza here.

Unirep & TAPE – highlight perfmrance on GFP prediction (mutant representation) 

## Approach

PPIDB:

3-Way split vs. random sampling:

Models:





## Results



## Discussion

## Citations

## Notes:
*DB:* To replicated the PPIDB, HuRI psi files for H-I-05, H-II-14, and HuRI must be [downloaded](http://www.interactome-atlas.org/download) and placed in the HuRI folder in ppiDB.  Similarly, the HIPPIE dataset must be downloaded and placed in the HIPPIE folder.  These files were not uploaded due to size restraints. 

*Embeddings:* TAPE was used to create the unirep embeddings.  To replicate the CD-HIT cluster dataset embeddings exactly, install [tape](https://github.com/songlab-cal/tape) and use:
'tape-embed unirep ppiDB.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10 ' and 'tape-embed unirep fragoza.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10'

*CD-HIT:* Install CD-HIT from [cdhit](https://github.com/weizhongli/cdhit) and use the default 


