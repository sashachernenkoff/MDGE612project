# MDGE612 Project
Sasha Chernenkoff

Predicting phenotype from genotype: a comparison of lasso and ridge regularized regression

## Introduction
This study aims to model, evaluate, and compare the power of various machine
learning models to predict a phenotype given the genotype. The models to be
compared include multiple linear regression (without regularization), L1 (LASSO) regularization
and L2 (Ridge) regularization. JawaMix5 (Long et al., 2013) was used to adjust for population
structure and perform association mapping to identify candidate SNPs from a list of
214,553 SNPs.


## Methods

![](https://github.com/sashachernenkoff/MDGE612project/blob/main/img/MDGE612%20pipeline.png?raw=true)

### Quality control
Python 3.11 was used with packages pandas and numpy for preprocessing the input genotype 
and phenotype files for modelling in Jawamix5. Preprocessing included imputing any missing 
genotypes, recoding the genotypes, and removing individuals with no phenotype data.


### Genome-wide association mapping
Jawamix5 was used to determine p-values of all SNPs. This tool incorporates EMMAX, a method that 
utilizes a mixed linear model to adjust for population structure, as described by Kang et al. 
(2010). Notably, this model includes a process to estimate the kinship matrix's influence on the 
phenotypes.


### Regression models
Python 3.11 was used with packages pandas, numpy, and scikitlearn to perform the predictions. For 
regularized models, grid search was performed using GridSearchCV for regression coeffiecient 
tuning. 


### Annotation
To provide biological context to the SNPs and identify potential functional implications, all 
significant SNPs identified by Jawamix5, as well as the top 100 SNPs used in the regression models, 
were annotated using the genomic features defined in the GFF file.


## Results
We identified 1,254 candidate SNPs from genome-wide association mapping. The SNPs were then sorted 
in ascending order based on their p-value. The top 20, 40, 60, 80, and 100 SNPs were used to compare 
the performance of the three regression models (linear regression, L1 regularization, and L2 
regularization). The L1 regularized model provided the lowest MSE and highest R^2 value.

![](https://github.com/sashachernenkoff/MDGE612project/blob/main/img/model_eval.png?raw=true)


### Analysis of top SNPs with GFF model
The majority of significant SNPs identified by Jawamix5 were coding region variants, particularly 
protein coding variants. Other minor frequency features included snRNA, ncRNA, 5' UTR, and 3' UTR. 
A similar distribution of features was observed in the top 100 most significant SNPs.

![](https://github.com/sashachernenkoff/MDGE612project/blob/main/img/annos.png?raw=true)


## Supplementary information
The additional code has been attached as python and unix command files.


## References
Atwell, S. *et al*. Genome-wide association study of 107 phenotypes in *rabidopsis thaliana* 
inbred lines. *Nature 2010 465:7298* 465, 627–631 (2010). https://doi.org/10.1038/nature08800

Kang, H. M. *et al*. Variance component model to account for sample structure in genome-wide 
association studies. *Nat Genet* 42, 348–354 (2010). https://doi.org/10.1038/ng.548

Long, Q. *et al*. JAWAMix5: an out-of-core HDF5-based java implementation of whole-genome 
association studies using mixed models. *Bioinformatics* 29, 1220–1222 (2013). 
https://doi.org/10.1093/bioinformatics/btt122