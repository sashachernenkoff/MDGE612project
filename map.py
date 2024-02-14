import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy.stats as stats
import statsmodels.api as sm

from utils import plot_manhattan, recode, ld_trim_2kb

# Significance level
sig = 0.05/214553

# Read in the genotypes and phenotypes, recode genotypes
# geno_raw = pd.read_csv('call_method_54.tair9.FT10.csv')
# geno = recode(geno_raw)
# geno.to_csv('coded_call_method_54.tair9.FT10.csv',index=False)
geno = pd.read_csv('coded_call_method_54.tair9.FT10.csv')
pheno = pd.read_csv('FT10.txt', sep='\t')

# Identify individuals with both genotype and phenotype data
geno.columns = geno.columns.astype(str)
com = pd.concat([geno, pd.DataFrame(np.nan, index=[len(geno)], columns=geno.columns)])

# Update the phenotype data in com
for index, row in pheno.iterrows():
    id_str = str(int(row['ecotype_id']))  # Why was this converting to a float?
    if id_str in com.columns:
        com.at[len(geno), id_str] = row[pheno.columns[1]]
    # else:
    #     print(f"ID {id_str} not found in com columns")

# Remove all individuals for which there is no phenotype data in com and separate out
# Transpose genotype dataframe to perform association mapping
com_fil = com.loc[:, com.loc[com.index[-1]].notna()]
geno_fil = com_fil.iloc[:-1].T
pheno_fil = com_fil.iloc[-1:].reset_index(drop=True)

# Calculate p-value for each snp
ps = []
y = pd.to_numeric(pheno_fil.iloc[0, :], errors='coerce')

# In pandas, using double square brackets around a column name returns a
# DataFrame with just that column, rather than a Series which would be
# returned by using single square brackets
for i, snp in enumerate(geno_fil.columns):
    X = geno_fil[[snp]]
    X = sm.add_constant(X)
    model = sm.OLS(y, X).fit()
    ps.append(model.pvalues[snp])
    if i % 10000 == 0:
        print(f'    processing snp {i}')

ps_log10 = [-np.log10(x) for x in ps]

results = pd.DataFrame({
    'Chromosome': geno['Chromosome'],
    'Positions': geno['Positions'],
    'p': ps,
    '-log10(p)': ps_log10
})

results = results.dropna()
results.to_csv('results.csv')
plot_manhattan(results, sig)

ssnps = results[results['p'] < sig]

results = pd.read_csv('results.csv',index_col=0)

# Correct for LD
snps_trim = ld_trim_2kb(results)
ssnps_trim = snps_trim[snps_trim['p'] < sig]
ssnps_trim.to_csv('ssnps.csv')

# Calculate the genomic inflation factor
# (estimates the amount of inflation by comparing
# observed test statistics across all genetic variants
# to those expected under the hypothesis of no effect
observed_chi2 = -2 * np.log(results['p'])
expected_median_chi2 = stats.chi2.ppf(0.5, df=1)
gif = np.median(observed_chi2) / expected_median_chi2

# Plot p-value distribution (QQ plot)
stats.probplot(observed_chi2, dist="chi2", sparams=(1,), plot=plt)
plt.title('QQ Plot')
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Observed Chi-Square Test Statistics')
plt.show()

# Perform PCA
geno_T = geno.iloc[:,2:].T

# Standardize the data
# It's important for PCA since it's sensitive to variances of your variables
geno_stand = StandardScaler().fit_transform(geno_T)

# typically a small number of PCs is sufficient for population structure
# e.g., 10 to 20 components
n = 10

pca = PCA(n_components=n)
pcs = pca.fit_transform(geno_stand)
pca_df = pd.DataFrame(data=pcs, columns=['PC' + str(i) for i in range(1, n+1)])

print(pca.explained_variance_ratio_)
print(sum(pca.explained_variance_ratio_))

# Scree plot
plt.plot(range(1,n+1), pca.explained_variance_ratio_, color='black', lw=1, marker='o',
         markersize=4, markerfacecolor='white')
plt.xlabel('PC number')
plt.ylabel('Variance (%)') # for each component
plt.title('Scree Plot')
plt.show()

# Plot PC1 and PC2
plt.scatter(pca_df['PC1'], pca_df['PC2'], c='black', marker='.')
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.title('PCA')
plt.show()


# take the 100 most significant SNPs/ or lowest 1 or 0.1 % 0.5%
