# embeddingsPPI: Binary protein-protien interaction prediction with unsupervised learning embeddings

Protein-protein interactions (PPIs) are physical contacts between two or more proteins and are typically essential for carrying out normal protein functions in the cell. Variation in genetic coding sequences can result in missense mutations in protein amino acid sequence, potentially disrupting PPIs and hindering or destroying normal protein function.  As such, reliably predicting how a missense mutation will affect an interaction between two proteins could is an interesting deep learning task to both better understand protein behavior and to differentiate benign and pathogenic mutations.  There have been multiple high-performance deep learning models released for PPI prediction, but these models perform poorly when predicting mutant PPI effects.  For this project, I evaluated unsupervised protein embedding models ability to improve mutant PPI prediction peformance. 

## Introduction


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
`tape-embed unirep ppiDB.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10` and `tape-embed unirep fragoza.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10`

*CD-HIT:* Install CD-HIT from [cdhit](https://github.com/weizhongli/cdhit) and use the default settings with `-c 0.7`when running the clustering program using a fasta file containing all sequences between 50 adn 2000 AA long.


