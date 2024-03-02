import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.rcParams['font.family'] = 'serif'

# A function for creating a manhattan plot from the output of Jawamix5
def plot_manhattan(snps, sig):
    snps = snps.copy()
    chrom_colors = ['darksalmon', 'cadetblue', 'goldenrod', 'rebeccapurple', 'olive']
    color_map = {chrom: chrom_colors[(chrom - 1) % len(chrom_colors)] for chrom in snps['chr'].unique()}
    snps['color'] = snps['chr'].map(color_map)

    # Create the Manhattan plot
    plt.figure(figsize=(8, 3))
    plt.scatter(snps.index, snps['-log10(p)'], c=snps['color'], alpha=0.6, s=10)

    # Add a horizontal line for the significance threshold
    plt.axhline(y=-np.log10(sig), color='grey', linestyle='--')

    # Calculate the midpoints of each chromosome region for labeling
    # and add chromosome labels at their midpoints
    midpoints = np.repeat(0, 5)
    for i in range(0,5):
        chrom = snps[snps['chr'] == i + 1]
        midpoints[i] = (chrom.index[0] + chrom.index[len(chrom)-1])/2

    # Add labels
    plt.title('Manhattan Plot')
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.xticks(midpoints,[1,2,3,4,5])
    plt.ylim(bottom=2.25)
    plt.xlim(left=0,right=len(snps))

    # Show the plot
    plt.savefig('data/jawamix5_manhattan_plot.pdf')
    plt.show()


# Read in data
colnames = ['chr', 'location', 'pvalue', 'AdjustedR2', 'coefficient', 'Sd_Err', 'MAF_count']
snps = pd.read_csv('data/emmax_out/EMMAX.0_5_FT10.top', skiprows=2, names=colnames)

# Calculate -log10(p)
snps['-log10(p)'] = [-np.log10(x) for x in snps['pvalue']]

# Plot manhattan
plot_manhattan(snps, 0.05)
