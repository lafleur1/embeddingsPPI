# embeddingsPPI: Binary protein-protien interaction prediction with unsupervised learning embeddings

Protein-protein interactions (PPIs) are physical contacts between two or more proteins and are typically essential for carrying out normal protein functions in the cell. Variation in genetic coding sequences can result in missense mutations in protein amino acid sequence, potentially disrupting PPIs and hindering or destroying normal protein function.  As such, reliably predicting how a missense mutation will affect an interaction between two proteins could is an interesting deep learning task to both better understand protein behavior and to differentiate benign and pathogenic mutations.  There have been multiple high-accuracy deep learning models released for PPI prediction, but these models performed poorly when predicting mutant PPI effects.  For this project, I evaluated unsupervised protein embedding models ability to improve mutant PPI prediction peformance. 

## Introduction
Variation in coding sequence (CDS) regions can greatly impact protein function and has the potential to disrupt PPIs. Recent work by (Fragoza et al, 2019) using the ExAC dataset of CDS variants from over 60K human exomes has estimated that 10.5% of an individual’s missense mutations disrupt PPIs, and indicated that most PPI disrupting variants do not affect overall protein stability but instead affect specific PPI interfaces by local structural disruptions.  Additionally, it has been shown that two-thirds of known pathogenic variants disrupt PPIs, with most of these disrupting only subsets of normal protein interaction profiles (Navío et al., 2019). Learning to reliably predict how small primary sequence changes can disrupt specific PPIs thus becomes an important step towards identifying putative pathogenic variants, and understanding how wild type protein interaction profiles can be disrupted through design efforts.

Several groups have developed predictive models that integrate protein structure information with information about measured protein-protein interactions (PPIs) to assess the likely impact of a variant (Engin, Kreisberg and Carter, 2016; Zhou et al., 2020), and multiple sequence only PPI binary predictors exist (Z. Chen et al., 2019). However, these models still have limited predictive power on protein families not seen during training, mainly due to limited training data and difficulties learning the complex sequence-structure relationships necessary for interaction predictions.

## Related Work

The best performing binary PPI prediction model to date is PIPR from (Chen et al., 2020), a describe here.  Describe embedding used.  Builds on similar ideas from (X), using a siamese CNN network to encode both one-hot encoded protein sequences
PIPR.  Another interesting PPI model is (NAME HEre), (MuPIPR), which 

Introduce Fragoza dataset. Show poor perofrmance on Fragoza here.

Unirep & TAPE – highlight perfmrance on GFP prediction (mutant representation) 

## Approach

### PPIDB
A SQL database of non-redundant PPIs was constructed by combining 

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


