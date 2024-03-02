import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.family'] = 'serif'

colnames = ['Chromosome', 'Positions', 'pvalue', 'AdjustedR2', 'coefficient', 'Sd_Err', 'MAF_count']
all_snps = pd.read_csv('data/emmax_out/EMMAX.0_5_FT10.top', skiprows=2, names=colnames)
top_snps = pd.read_csv('data/ssnps_jawamix5.csv')

gff_columns = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']
gff = pd.read_csv('gene_model.gff', sep='\t', names=gff_columns, comment='#')


def find_types(chromosome, position):
    chr = f'Chr{int(chromosome)}'
    chr_feat = gff[gff['Chromosome'] == chr]
    feat = chr_feat[(chr_feat['Start'] <= position) & (chr_feat['End'] >= position)]
    return ";".join(feat['Type']) if not feat.empty else float('NaN')


# Annotate all snps
all_snps['FeatureTypes'] = all_snps.apply(lambda row: find_types(row['Chromosome'], row['Positions']), axis=1)
all_snps.to_csv('data/all_annotated_ssnps.csv')

# Annotate the top snps
top_snps['FeatureTypes'] = top_snps.apply(lambda row: find_types(row['Chromosome'], row['Positions']), axis=1)
top_snps.to_csv('data/top_annotated_ssnps.csv')


def count_annotations(df):
    all_types = ';'.join(df['FeatureTypes'].dropna()).split(';')
    all_types = [atype for atype in all_types if atype.lower() != "chromosome"]

    type_counts = pd.Series(all_types).value_counts()
    return type_counts


# Count annotations for snp set
all_snps_counts = count_annotations(all_snps)
top_snps_counts = count_annotations(top_snps)

# Plot annotation counts
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

# All SNPs
all_snps_counts.plot(kind='barh', ax=ax[0], color='darksalmon')
ax[0].set_title('Annotations for all SNPs')
ax[0].set_xlabel('Frequency')
ax[0].set_ylabel('Type of annotation')

# Top SNPs
top_snps_counts.plot(kind='barh', ax=ax[1], color='cadetblue')
ax[1].set_title('Annotations for top 100 SNPs')
ax[1].set_xlabel('Frequency')
# ax[1].set_ylabel('Type of annotation')  # Already set by the first plot

plt.tight_layout()
plt.savefig('img/annos.pdf')
plt.show()