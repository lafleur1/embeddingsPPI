# embeddingsPPI: Binary protein-protien interaction prediction with unsupervised learning embeddings

Protein-protein interactions (PPIs) are physical contacts between two or more proteins and are typically essential for carrying out normal protein functions in the cell. Variation in genetic coding sequences can result in missense mutations in protein amino acid sequence, potentially disrupting PPIs and hindering or destroying normal protein function.  As such, reliably predicting how a missense mutation will affect an interaction between two proteins could is an interesting deep learning task to both better understand protein behavior and to differentiate benign and pathogenic mutations.  There have been multiple high-accuracy deep learning models released for PPI prediction, but these models performed poorly when predicting mutant PPI effects.  For this project, I attempted to evaluate unsupervised protein embedding models from the Tasks Assessing Protein Embedding (TAPE) benchnark to improve mutant PPI prediction peformance. 

## Introduction
Variation in coding sequence (CDS) regions can greatly impact protein function and has the potential to disrupt PPIs. Recent work by (Fragoza et al, 2019) using the ExAC dataset of CDS variants from over 60K human exomes has estimated that 10.5% of an individual’s missense mutations disrupt PPIs, and indicated that most PPI disrupting variants do not affect overall protein stability but instead affect specific PPI interfaces by local structural disruptions.  Additionally, it has been shown that two-thirds of known pathogenic variants disrupt PPIs, with most of these disrupting only subsets of normal protein interaction profiles (Navío et al., 2019). Learning to reliably predict how small protein sequence changes can disrupt specific PPIs thus becomes an important step towards identifying putative pathogenic variants, and understanding how wild type protein interaction profiles can be disrupted through protein design efforts.

Several groups have developed predictive models that integrate protein structure information with information about measured protein-protein interactions (PPIs) to assess the likely impact of a variant (Engin, Kreisberg and Carter, 2016; Zhou et al., 2020), and multiple sequence only PPI binary predictors exist (Z. Chen et al., 2019). However, these models still have limited predictive power on protein families not seen during training, mainly due to limited training data and difficulties learning the complex sequence-structure relationships necessary for interaction predictions.

## Related Work
The prediction of PPI occurrence is a widely attempted deep learning task, with recent high accuracy and precision models adapting natural language processing (NLP) methods for protein sequence representations (M. Chen et al., 2019; Nambiar et al., 2020).  Networks have also been trained to predict interaction type, binding affinity, and interaction sites from hand-picked features, sequence, and crystal structures (Geng et al., 2019; M. Chen et al., 2019; Gainza et al., 2020).

The best performing binary PPI prediction model to date is PIPR from (Chen et al., 2020), a deep residual recurrent siamese CNN.  This network embeds protein sequences using amino-acid co-occurence similairty and a one-hot encoded amino acid sequence.  The embedded sequences are fed into the siamese residual RCNN portion wtih shared parameters to create a sequence embedding vector for each potentially interacting sequence.  These vectors are combined with element-wise multiplication, before being fed through a linear netowrk portion for interaction prediction.  Another state-of-the-art PPI property predictor model is , (MuPIPR), which predicts binding affinity change between a wild type protein pair and a mutant protien pair and uses a pre-trained bidirectional language model (BiLM) to embed all four protein input sequences into the siamese netowrk which is then fed into a residuan CNN.  Therefore, the use of NLP models for PPI prediction is well-established.

A dataset which can used to judge the effectiveness of these state-ofthe-art models for predicting how protein sequence changes can disrupt specific PPIs was recently collected by (Fragoza et al., 2019).  Interaction profiles for 4797 SNV-interaction pairs were measured using yeast two-hybrid assays, a common PPI measurement method using modified yeast to detect if two proteins interact.  Of these, a subset can be used to evaluate PIPR and MuPIPR's ability to predict mutatnt PPIs a shown in **Figure 1**. As it can be seen, these models fail to separate mutatiosn whic 

Unirep & TAPE – highlight perfmrance on GFP prediction (mutant representation) 

## Approach

### PPIDB
A SQL database of non-redundant PPIs was constructed by combining 

3-Way split vs. random sampling:

Models:





## Results
As the main goal of this project was to determine if further pursuit of UniRep and Bert embeddngs should be pursued to improve mutant 


## Discussion
While the 

## Citations

## Notes:
*DB:* To replicated the PPIDB, HuRI psi files for H-I-05, H-II-14, and HuRI must be [downloaded](http://www.interactome-atlas.org/download) and placed in the HuRI folder in ppiDB.  Similarly, the HIPPIE dataset must be downloaded and placed in the HIPPIE folder.  These files were not uploaded due to size restraints. 

*Embeddings:* TAPE was used to create the unirep embeddings.  To replicate the CD-HIT cluster dataset embeddings exactly, install [tape](https://github.com/songlab-cal/tape) and use:
`tape-embed unirep ppiDB.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10` and `tape-embed unirep fragoza.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10`

*CD-HIT:* Install CD-HIT from [cdhit](https://github.com/weizhongli/cdhit) and use the default settings with `-c 0.7`when running the clustering program using a fasta file containing all sequences between 50 adn 2000 AA long.


