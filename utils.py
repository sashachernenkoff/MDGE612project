import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats


# A function for reading in a genotype file and recoding the genotype to
# binary data, where 0 represents the reference allele and 1 represents
# an alternative allele.
def recode(genotypes):
    geno_coded = genotypes
    ref_alleles = pd.DataFrame({'allele': ['N'] * len(genotypes)}, index=genotypes.index)

    # Find the most common allele
    for i, row in geno_coded.iterrows():
        ref_alleles.iloc[i] = row.mode().astype(str)

    for i,row in geno_coded.iterrows():
        for j in range(2, len(row)):
            if row[j] == ref_alleles.iloc[i].item():
                geno_coded.at[row.name, geno_coded.columns[j]] = 0
            else:
                geno_coded.at[row.name, geno_coded.columns[j]] = 1
        if i % 10000 == 0:
            print(f'    recoding row {i}')

    return geno_coded


# A function for creating a manhattan plot from a data frame containing snps
# pre-sorted by 'Chromosome', 'Position', with pre-calculated '-log(10)p'
def plot_manhattan(snps, sig):
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
    plt.savefig('plot.pdf')
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
