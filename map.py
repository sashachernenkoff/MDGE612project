import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm


# A function for creating a manhattan plot from a data frame containing snps
# pre-sorted by 'Chromosome', 'Position', with pre-calculated '-log(10)p'
def plot_manhattan(snps, sig):
    snps = snps.copy()
    chrom_colors = ['darksalmon', 'cadetblue', 'goldenrod', 'rebeccapurple', 'olive']
    color_map = {chrom: chrom_colors[(chrom - 1) % len(chrom_colors)] for chrom in snps['Chromosome'].unique()}
    snps['color'] = snps['Chromosome'].map(color_map)

    # Create the Manhattan plot
    plt.figure(figsize=(12, 6))
    plt.scatter(snps.index, snps['-log10(p)'], c=snps['color'], alpha=0.6, s=10)

    # Add a horizontal line for the significance threshold
    plt.axhline(y=-np.log10(sig), color='grey', linestyle='--')

    # Calculate the midpoints of each chromosome region for labeling
    # and add chromosome labels at their midpoints
    midpoints = np.repeat(0, 5)
    for i in range(0,5):
        chrom = snps[snps['Chromosome'] == i + 1]
        midpoints[i] = (chrom.index[0] + chrom.index[len(chrom)-1])/2

    # Add labels
    plt.title('Manhattan Plot')
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.xticks(midpoints,[1,2,3,4,5])

    # Show the plot
    plt.savefig('data/manhattan_plot.pdf')
    plt.show()


# A function that selects the snp with the most significant p-value for every
# 2kbp bin in a list of snps
def ld_trim_2kb(snps):
    snps_trim = pd.DataFrame(columns=snps.columns)

    for i in range(1, 6):
        print(f'chrom {i}')
        chrom = snps[snps['Chromosome'] == i]
        for j in range(2000, chrom['Positions'].max(), 2000):
            if j % 10000000 == 0:
                print(f'    chrom {i}, pos {j}')
            bin = chrom[(chrom['Positions'] >= j - 2000) & (chrom['Positions'] < j)]
            if not bin.empty:
                # .iloc[0] to take just one snp (if more than one)
                snps_trim.loc[len(snps_trim.index)] = bin[bin['p'] == bin['p'].min()].iloc[0]
        # take the snp from the final bin here!

    snps_trim['Chromosome'] = snps_trim['Chromosome'].astype(int)
    return snps_trim


# Significance level (Bonferroni correction)
sig = 0.05/214553

# Read in the genotypes and phenotypes, recode genotypes
geno = pd.read_csv('data/filtered.coded_call_method_54.tair9.FT10.csv')
pheno = pd.read_csv('data/filtered.FT10.txt')

# Calculate p-value for each snp
ps = []
y = pd.to_numeric(pheno['5_FT10'], errors='coerce')

# In pandas, using double square brackets around a column name returns a
# DataFrame with just that column, rather than a Series which would be
# returned by using single square brackets
X_geno = geno.drop(['Chromosome', 'Positions'], axis=1).T
for i, snp in enumerate(X_geno.columns):
    X = X_geno[[snp]].reset_index(drop=True)
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

# Drop snps which don't have a p-value
results = results.dropna() # (214553 x 4)
results.to_csv('data/map_results.csv') # (214219 x 4)
# results = pd.read_csv('data/map_results.csv', index_col=0)

plot_manhattan(results, sig)

# Correct for LD
ssnps = results[results['p'] < sig] # 4180 significant snps before trim
snps_trim = ld_trim_2kb(results)
ssnps_trim = snps_trim[snps_trim['p'] < sig]
ssnps_trim.to_csv('data/ssnps.csv') # 2764 significant snps after trim
# ssnps = pd.read_csv('data/ssnps.csv', index_col=0)

# Calculate the genomic inflation factor
# (estimates the amount of inflation by comparing
# observed test statistics across all genetic variants
# to those expected under the hypothesis of no effect
observed_chi2 = -2 * np.log(results['p'])
expected_median_chi2 = stats.chi2.ppf(0.5, df=1)
gif = np.median(observed_chi2) / expected_median_chi2
print(f'genomic inflation factor: {gif}') # gif = 9.111

# Plot p-value distribution (QQ plot)
stats.probplot(observed_chi2, dist="chi2", sparams=(1,), plot=plt)
plt.title('QQ Plot')
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Observed Chi-Square Test Statistics')
plt.show()

# Need to handle population stratification