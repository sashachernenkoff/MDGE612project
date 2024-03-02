import numpy as np
import pandas as pd


# A function for reading in a genotype file and recoding the genotype to
# binary data, where 0 represents the reference (most common) allele and 1
# represents an alternative allele.
def recode(genotypes):
    genotype_data = genotypes.iloc[:, 2:]
    refs = genotype_data.mode(axis=1)[0]
    matches_ref = genotype_data.apply(lambda x: x == refs.loc[x.index], axis=0)
    matches_ref = ~matches_ref
    recoded = pd.concat([genotypes.iloc[:, :2], matches_ref.astype(int)], axis=1)
    return recoded


# Read in the genotypes and phenotypes, recode genotypes
geno_raw = pd.read_csv('data/call_method_54.tair9.FT10.csv')

# Check for and remove NAs
# missing = ~geno_raw.iloc[:,2:].isin(['A', 'C', 'T', 'G'])
# print(missing.sum().sum())
# No missing genotypes

geno = recode(geno_raw)
geno.to_csv('out/coded_call_method_54.tair9.FT10.csv',index=False)
# geno = pd.read_csv('data/coded_call_method_54.tair9.FT10.csv')
pheno = pd.read_csv('data/FT10.txt', sep='\t')

# Identify individuals with both genotype and phenotype data
com = pd.concat([geno.iloc[:,2:], pd.DataFrame(np.nan, index=[len(geno)], columns=geno.columns)])

for i, row in pheno.iterrows():
    sid = str(int(row['ecotype_id']))
    if sid in com.columns:
        com.at[len(geno), sid] = row[pheno.columns[1]]

# Remove all individuals for which there is no phenotype data
com = com.dropna(axis=1)

# Sort and save genotype and phenotype files
com = com.reindex(sorted(com.columns, key=lambda x: int(x)), axis=1)
geno_fil = pd.concat([geno_raw[['Chromosome', 'Positions']], com.drop(len(com)-1).astype(int)], axis=1)
pheno_fil = pd.DataFrame({
    'ecotype_id':com.columns,
    '5_FT10': com.iloc[-1]}).reset_index(drop=True)

geno_fil.to_csv('out/filtered.coded_call_method_54.tair9.FT10.csv',index=False)
pheno_fil.to_csv('out/filtered.FT10.txt',index=False)