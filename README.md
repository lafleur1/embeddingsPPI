# embeddingsPPI: Binary protein-protien interaction prediction with semi-supervised embeddings

Protein-protein interactions (PPIs) are physical contacts between two or more proteins and are typically essential for carrying out normal protein functions in the cell. Variation in genetic coding sequences can result in missense mutations in protein amino acid sequence, potentially disrupting PPIs and hindering or destroying normal protein function.  As such, reliably predicting how a missense mutation will affect an interaction between two proteins could is an interesting deep learning task to both better understand protein behavior and to differentiate benign and pathogenic mutations.  There have been multiple high-accuracy deep learning models released for PPI prediction, but these models performed poorly when predicting mutant PPI effects.  For this project, I attempted to evaluate unsupervised protein embedding models from the Tasks Assessing Protein Embedding (TAPE) benchnark (Rao et al., 2019) for binary PPI predictions, and to create a large high-quality dataset for use in furhter PPI modeling attempts to improve mutant PPI prediction. 

## Introduction
Variation in coding sequence (CDS) regions can greatly impact protein function and has the potential to disrupt PPIs. Recent work by (Fragoza et al, 2019) using the ExAC dataset of CDS variants from over 60K human exomes has estimated that 10.5% of an individual’s missense mutations disrupt PPIs, and indicated that most PPI disrupting variants do not affect overall protein stability but instead affect specific PPI interfaces by local structural disruptions.  Additionally, it has been shown that two-thirds of known pathogenic variants disrupt PPIs, with most of these disrupting only subsets of normal protein interaction profiles (Navío et al., 2019). Learning to reliably predict how small protein sequence changes can disrupt specific PPIs thus becomes an important step towards identifying putative pathogenic variants, and understanding how wild type protein interaction profiles can be disrupted through protein design efforts.

Several groups have developed predictive models that integrate protein structure information with information about measured protein-protein interactions (PPIs) to assess the likely impact of a variant (Engin, Kreisberg and Carter, 2016; Zhou et al., 2020), and multiple sequence only PPI binary predictors exist (Z. Chen et al., 2019). However, these models still have limited predictive power on protein families not seen during training, mainly due to limited training data and difficulties learning the complex sequence-structure relationships necessary for interaction predictions.

## Related Work
The prediction of PPI occurrence is a widely attempted deep learning task, with recent high accuracy and precision models adapting natural language processing (NLP) methods for protein sequence representations (M. Chen et al., 2019; Nambiar et al., 2020).  Networks have also been trained to predict interaction type, binding affinity, and interaction sites from hand-picked features, sequence, and crystal structures (Geng et al., 2019; M. Chen et al., 2019; Gainza et al., 2020).

The best performing binary PPI prediction model to date is PIPR from (Chen et al., 2020), a deep residual recurrent siamese CNN.  This network embeds protein sequences using amino-acid co-occurence similairty and a one-hot encoded amino acid sequence.  The embedded sequences are fed into the siamese residual RCNN portion wtih shared parameters to create a sequence embedding vector for each potentially interacting sequence.  These vectors are combined with element-wise multiplication, before being fed through a linear netowrk portion for interaction prediction.  Another state-of-the-art PPI property predictor model is MuPIPR, which predicts binding affinity change between a wild type protein pair and a mutant protien pair and uses a pre-trained bidirectional language model (BiLM) to embed all four protein input sequences into the siamese netowrk which is then fed into a residuan CNN (Zhou et al., 2020).  Therefore, the use of NLP models for PPI prediction is well-established.

A dataset which can used to judge the effectiveness of these state-ofthe-art models for predicting how protein sequence changes can disrupt specific PPIs was recently collected by (Fragoza et al., 2019).  Interaction profiles for 4797 SNV-interaction pairs were measured using yeast two-hybrid assays, a common PPI measurement method using modified yeast to detect if two proteins interact.  Of these, a subset can be used to evaluate PIPR and MuPIPR's ability to predict mutatnt PPIs a shown in **Figure 1**. As it can be seen, these models fail to separate PPI disrutpive mutants from non-disruptive mutants in the Fragoza dataset. 

FIGURE 1 GOES EHRE 

Recently, an attempt was made to standardize and evaluate the effecivness of unsupervised and semi-supervised learning techniques on several protien biology deep learning tasks, Tasks Assessing Protein Embeddings (TAPE) (Rao et al., 2019).  TAPE foudn that semi-supervised NLP pretraining methods helped improve performance on most of the tasks tested, and assessed a transformer, LSTM, ResNet, and a multiplcative LSTM (mLSTM) for structure, evolutionary, and engieerng prediciton tasks. However, there was no best performing pretrained NLP method across all assessed tasks.  Two of the semi-supervised protein sequence models were made publically avaible form TAPE, bert_base, a transformer, and UniRep, a mLSTM.  As these are different than the NLP methods employeed in PIPR and MuPIPR, tehy could potentially improve PPI prediciton abilities (Chen et al., 2020, Zhou et al., 2020).

## Approach

### Evaluating current models on the Fragoza dataset
PIPR was cloned from [seq_ppi](https://github.com/muhaochen/seq_ppi), and all Fragoza datset interactions where both partners had a maximum length of <= 600 amino acids were evaluted.  Code to reproduce the PIPR portion of **Figure 1** is located in *currentModelFragozaComparisons/PIPR/*.  MuPIPR was cloned from [MuPIPR](https://github.com/guangyu-zhou/MuPIPR), and again all Fragoza datset interactions where both partners had a maximum length of <= 600 amino acids were evaluted.  Code to reproduce the MuPIPR portion of **Figure 1** is located in *currentModelFragozaComparisons/MuPIPR/*. (Interactions where both partners had a maximum length of <= 600 amino acids were used in this evaluation as it was the maximum input length possible between the two models). 

### ppiDB
A SQL database of non-redundant PPIs, ppiDB, was constructed by combining several human PPI databases and high-throughput datasets. A simplified schema of the database is shown in **Figure 2**.  Datasets were downloaded from their respective websites (given in links below) and gene IDs, transcript IDs, or Uniprot IDs were retrieved for available sequences with the Uniprot *Retrieve/ID mapping tool*.  Subsets of ENSEMBL sequences were downloaded as needed if a specific transcript ID was provided in the database.  Deleted and deprecated sequences from Uniprot were not retrieved to be incorporated into ppiDB.  The databases and datasets incorporated were:

- The Human Reference Interactome (HuRI): (Luck et al., 2020) A high-throughput human PPI dataset, part of an ongoing attempt to generate a complete draft of the human interactome.  Used as a source of positive interactions (those in the dataset) and negative itneractions (interactions which were screened in the all-by-all Y2H experiments, but an interaction was not detected).
- Human Integrated Protein Protein Interaction Reference (HIPPIE) v2.2: (Alanis-Lobato et al., 2017) A confidence scored and functionally annotated human PPIs database. Used as a source of positive interactions.
- HitPredict v4: (Lopez et al., 2015) A annotated and scored database of PPIs with sources from several data repositories. Used as a source of positive interacitons
- Negatome 2.0: (Blohm et al., 2014) A database of protein pairs deemed unlikeky to interact, either curated from the literature or determined from protein structures. Used as source of negative interactions.

FIGURE 2 GOES EHRE 

### Dataset extraction from ppiDB
TAPE embeddings were trained on sequences with a maximum lenght of 2000 amino acids.  All sequences with lenght <= 2000 were extracted from ppiDB and placed in a fasta file.  They were then clustered at a 70% sequence identity cutofff with CD-HIT (Limin et al., 2012), a common technique in biology deep learning models to prevent leak between training and evaluation datasets due to the high level of similarity between protein sequences.  This resulted in 16614 representative cluster sequences.  Positive interactions from HuRI, HIPPIE, and HitPredict were gathered from ppiDB between cluster representatives, as well as Negatome negatives.  This resulted in 279,775 positive interactions and 918 negative interactions between the cluster sequences.  An equal number of negative interactions were then added to ppiDB from non-interacting pairs in HuRI to create a balanced dataset for network training, resulting in a total dataset of 279,775 positive interactions and 279,775 negative interactions.

### TAPE embeddings
Following cluster sequence selection, pre-trained unsupervised models UniRep and bert-base from TAPE were used to process cluster sequences. To replicate the CD-HIT cluster dataset embeddings exactly, install [tape](https://github.com/songlab-cal/tape) and use:
`tape-embed unirep ppiDB.fasta ppiDB.npz babbler-1900 --tokenizer unirep --batch_size = 10`.  Bert-base was processed in `runBert.ipynb` above.

### 3-Way split vs. random sampling:
As the most interesting case in PPI prediction is how well the model performs on predicting pairs which are highly dissimilar from the training set, two different versions of a train/validation/test split were conducted.  In one, a fraction of dataset proteins were chosen to be part of a held out set from training from the cluster proteins, resulting in interactions were both proteins were seen during training, only one was seen during training and one was unseen, and pairs where both proteins were unseen during training.  In the other, training, validation, adn test set pairs were uniform randomly selected from the dataset. 

### Models:
As the main goal of this project was to create the ppiDB for future embedding attempts, and to investigate the initial effectiveness of the pre-trained TAPE embeding methods, minimal network hyperparameter exploration was completed and very simple linear networks were trained for a small number of iterations to determine if any loss decreases would occur in the intial epochs.  Both representation methods tested used the same simple netowrk structure of the TAPE representations of interaction proteins fed into a 300 neuron fully connected layer followed by a 100 neuron fully connected layer iwth a Relu activation, followed by a single neuron with binary cross-entropy loss.  All bert_base models were trained for 5 epochs with a batch size of 32, using the Adam optimizer with a learning rate of 0.001.  The UniRep mdoels were were trained for 3 epochs with a batch size of 32, using the Adam optimizer with a learning rate of 0.001.  

Three methods of combinign the interaction protein partners were used: concatenation of the two input sequence representations, element-wise multiplicaiton of the input reseprentations, and differnce fo the input representations. 

## Results
The created database, ppiDB, is nearly an order of magnitude larger than the datasets used in training models such as PIPR and MuPIPR (Chen et al., 2020, Zhou et al., 2020), and can be a resource moving forward with testing new protein sequence representation methods to improve PPI prediction.  It contains 22561 protein sequences participating in 477,278 experimentally supported human PPIs.  This is a stark increase compared to the benchmark datasets typically used in the literature for binary PPI prediction, with the largest PPI training dataset used in PIPR containing only 11,529 proteins parcipating in 32,959 positive interactions and 32,959 negative interactions.  

All six model structures showed close to random performance in accuracy, as can be seen in Table 1. 

Table 1. Bert-base and UniRep Model Performance
| Attempt | #1  | #2  | #3 | #4 |
| :-----: | :-: | :-: | :-: | :-:|
| Seconds | 301 | 283 | x | y|
| Seconds | 301 | 283 | x | y|
| Seconds | 301 | 283 | x | y|

## Discussion
As the main goal of this project was to create a database that could be used in 
determine if further pursuit of UniRep and Bert embeddngs could improve mutant PPI prediction, it was useful in that it showed that the pooled adn average otuputs for these NLP based models are not appropriate for PPI prediction tasks.  Or, at the very least, they require a more specialized approach than was attemtped here.  

## Citations



## Notes:
*DB:* To replicated the PPIDB, HuRI psi files for H-I-05, H-II-14, and HuRI must be [downloaded](http://www.interactome-atlas.org/download) and placed in the HuRI folder in ppiDB.  Similarly, the HIPPIE dataset must be downloaded and placed in the HIPPIE folder.  These files were not uploaded due to size restraints. 

*CD-HIT:* Install CD-HIT from [cdhit](https://github.com/weizhongli/cdhit) and use the default settings with `-c 0.7`when running the clustering program using a fasta file containing all sequences between 50 adn 2000 AA long.

