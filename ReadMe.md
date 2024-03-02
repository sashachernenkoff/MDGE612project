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
and phenotype files for modelling in Jawamix5. Preprocessing included recoding the genotypes, 
removing individuals with no phenotype data, and imputing any missing genotypes.


### Genome-wide association mapping
Jawamix5 was used to obtain the significant SNPs and their p-values. This included converting 
the genotype file to hdf5, creating a kinship matrix, and performing the association mapping
(with EMMAX).


### Regression models
Python 3.11 was used with packages pandas, numpy, and scikitlearn to perform the predictions. For 
regularized models, grid search was performed using GridSearchCV for regression coeffiecient 
tuning. 


### Annotation
Describe annotation process


## Results
Jawamix5 was used to determine p-values of all SNPs. This tool incorporates EMMAX, a method that 
utilizes a mixed linear model to adjust for population structure, as described by Kang et al. 
(2010). Notably, this model includes a process to estimate the kinship matrix's influence on the 
phenotypes, a feature not present in other methods. We identified 1,254 candidate SNPs from genome-
wide association mapping. The SNPs were then organized in ascending order based on their minimum 
p-value. The top 20-100 SNPs were used to compare the performance of the three regression models.

![](https://github.com/sashachernenkoff/MDGE612project/blob/main/img/model_eval.png?raw=true)


### Analysis of top SNPs with GFF model
Insert discussion here

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