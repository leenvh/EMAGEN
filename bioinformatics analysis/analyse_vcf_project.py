
import numpy as np
import allel
import csv
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from matplotlib.patches import Ellipse
import pickle
from scipy.spatial.distance import squareform
from typing import List
import plotly.express as px
import pysam
import argparse
from scipy.stats import zscore
import matplotlib.pyplot as plt
import re
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from adjustText import adjust_text
from matplotlib.colors import LinearSegmentedColormap

#Prior to this script: 
#Filter out indels: bcftools view -v snps yourvcf.gz -Oz -o yourvcf.snps.vcf.gz


parser = argparse.ArgumentParser(description='Analyse VCF file')
parser.add_argument('vcf', type=str, help='VCF file')

args = parser.parse_args()

callset = allel.read_vcf(
    args.vcf, 
    fields=['samples', 'calldata/GT', 'calldata/DP', 'calldata/AD',
            'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
            'variants/QUAL', 'variants/ID', 'variants/FILTER_PASS']
)

#Depth filtering - set low depth genotypes to missing
if 'calldata/DP' in callset:
    dp = callset['calldata/DP']
    min_depth = 5
    
    # Create mask for low depth calls 
    low_depth_mask = (dp < min_depth)
    
    print(f"Genotype calls with DP < {min_depth}: {low_depth_mask.sum()} ({low_depth_mask.sum() / low_depth_mask.size * 100:.2f}%)")
    
    # Set low depth genotypes to missing (-1)
    gt_array = callset['calldata/GT'].copy()
    gt_array[low_depth_mask] = -1
    callset['calldata/GT'] = gt_array
else:
    print("Warning: No DP field found in VCF, skipping depth filtering")


def subset_callset(callset,variant_i,operation='extract'):
    if operation == 'exclude':
        mask = np.ones(len(callset['variants/POS']),dtype=bool)
        mask[variant_i] = False
    else:
        mask = np.zeros(len(callset['variants/POS']),dtype=bool)
        mask[variant_i] = True

    return {
        "samples":callset['samples'],
        "calldata/GT": callset['calldata/GT'][mask],
        "variants/ALT": callset['variants/ALT'][mask],
        "variants/CHROM": callset['variants/CHROM'][mask],
        "variants/FILTER_PASS": callset['variants/FILTER_PASS'][mask],
        "variants/ID": callset['variants/ID'][mask],
        "variants/POS": callset['variants/POS'][mask],
        "variants/QUAL": callset['variants/QUAL'][mask],
        "variants/REF": callset['variants/REF'][mask],
    }

# Select mitochondria only
#mito_idx = [i for i in range(len(callset['variants/CHROM'])) 
            #if callset['variants/CHROM'][i] in ['Pf_M76611']]
#print(f"Total variants before filtering: {len(callset['variants/POS'])}")
#print(f"Mitochondrial variants: {len(mito_idx)}")
#callset = subset_callset(callset, mito_idx, operation='extract')


vcf = pysam.VariantFile(args.vcf)
ann = {}
all_variant_types = [] 
for var in vcf:
    for d in var.info['ANN']:
        d = d.split('|')
        variant_type = d[1]
        all_variant_types.append(d[1])
        # Skip if not one of these types
        if variant_type not in ['missense_variant', 'synonymous_variant', 'upstream_gene_variant']:
            continue
            
        key = (var.chrom, var.pos, var.alts[0])
        ann[key] = (d[3].replace('gene:',''), d[10].replace('p.',''), variant_type)

from collections import Counter

print("\n=== All Variant Annotations (Before Any Filtering) ===")
variant_type_counts = Counter(all_variant_types)
print("\nVariant types and their frequencies:")
for variant_type, count in sorted(variant_type_counts.items(), key=lambda x: x[1], reverse=True):
    percentage = (count / len(all_variant_types)) * 100
    print(f"  {variant_type}: {count} ({percentage:.2f}%)")
print(f"\nTotal annotations: {len(all_variant_types)}")


variant_idx = [
    i for i in range(len(callset['variants/POS'])) 
    if (callset['variants/CHROM'][i], callset['variants/POS'][i], callset['variants/ALT'][i][0]) in ann
]
print("Number of variants: ",len(variant_idx))
callset = subset_callset(callset,variant_idx,operation='extract')

# add annotations to callset
callset['variants/GENE'] = []
callset['variants/AA'] = []
callset['variants/TYPE'] = []
for i in tqdm(range(len(callset['variants/POS']))):
    chrom = callset['variants/CHROM'][i]
    pos = callset['variants/POS'][i]
    alt = callset['variants/ALT'][i][0]
    callset['variants/GENE'].append(ann[(chrom,pos,alt)][0])
    callset['variants/AA'].append(ann[(chrom,pos,alt)][1])
    callset['variants/TYPE'].append(ann[(chrom,pos,alt)][2])




##### PCA 

#print(callset['variants/AA'])
g = allel.GenotypeArray(callset['calldata/GT'])
print(f"Starting dimensions: {g.shape[0]} variants x {g.shape[1]} samples")


# Calculate initial missingness
missing_per_variant = g.count_missing(axis=1) / g.shape[1]
missing_per_sample = g.count_missing(axis=0) / g.shape[0]

print(f"Variant missingness: min={missing_per_variant.min():.3f}, max={missing_per_variant.max():.3f}, mean={missing_per_variant.mean():.3f}")
print(f"Sample missingness: min={missing_per_sample.min():.3f}, max={missing_per_sample.max():.3f}, mean={missing_per_sample.mean():.3f}")

# Filter samples with high missingness
sample_missing_threshold = 0.5 
sample_mask = missing_per_sample < sample_missing_threshold
samples_removed = callset['samples'][~sample_mask]

print(f"Removing {len(samples_removed)} samples with >{sample_missing_threshold*100}% missing data:")
print(samples_removed)

# Apply sample filter
g_filtered = g.compress(sample_mask, axis=1)
samples_filtered = callset['samples'][sample_mask]

print(f"After sample filter: {g_filtered.shape[0]} variants x {g_filtered.shape[1]} samples")

# Filter variants with high missingness
variant_missing_threshold = 0.4  # Remove variants with >40% missing data
missing_per_variant_updated = g_filtered.count_missing(axis=1) / g_filtered.shape[1]
variant_mask_missing = missing_per_variant_updated < variant_missing_threshold

print(f"Removing {(~variant_mask_missing).sum()} variants with >{variant_missing_threshold*100}% missing data")

# Filter by Minor Allele Frequency (MAF)
g_temp = g_filtered.compress(variant_mask_missing, axis=0)
ac = g_temp.count_alleles()

# Calculate allele frequency
af = ac.to_frequencies()

# For biallelic sites, MAF is min(af[0], af[1])
# But some sites might be monomorphic after filtering, so handle that
maf = np.minimum(af[:, 0], af[:, 1]) if af.shape[1] > 1 else np.zeros(af.shape[0])

maf_threshold = 0.001 
variant_mask_maf = maf >= maf_threshold

print(f"Removing {(~variant_mask_maf).sum()} variants with MAF < {maf_threshold}")

# Combine variant filters and apply
variant_indices_missing = np.where(variant_mask_missing)[0]
variant_indices_final = variant_indices_missing[variant_mask_maf]

# Create final filtered genotype array
g_final = g_filtered.take(variant_indices_final, axis=0)

print(f"Final dimensions: {g_final.shape[0]} variants x {g_final.shape[1]} samples")

# Also filter the annotation data to match
callset_filtered = {
    'samples': samples_filtered,
    'variants/CHROM': callset['variants/CHROM'][variant_indices_final] if sample_mask.all() else np.array(callset['variants/CHROM'])[variant_indices_final],
    'variants/POS': callset['variants/POS'][variant_indices_final],
    'variants/REF': callset['variants/REF'][variant_indices_final],
    'variants/ALT': callset['variants/ALT'][variant_indices_final],
    'variants/GENE': [callset['variants/GENE'][i] for i in variant_indices_final],
    'variants/AA': [callset['variants/AA'][i] for i in variant_indices_final],
    'variants/TYPE': [callset['variants/TYPE'][i] for i in variant_indices_final],
}


# Convert to alternate allele counts and handle missing data
gn = g_final.to_n_alt(fill=-1)
print(f"Missing values in gn: {(gn == -1).sum()} ({(gn == -1).sum() / gn.size * 100:.2f}%)")


gn_imputed = gn.astype(float).copy()

for i in range(gn_imputed.shape[0]):
    variant_data = gn_imputed[i, :]
    missing_mask = variant_data == -1
    if missing_mask.any():
        non_missing = variant_data[~missing_mask]
        if len(non_missing) > 0:
            mean_val = non_missing.mean()
            gn_imputed[i, missing_mask] = mean_val
        else:
            # All missing - set to 0 
            gn_imputed[i, :] = 0

print(f"Missing values after imputation: {(gn_imputed == -1).sum()}")

np.savetxt(
    "gn_imputed_projects.csv",
    gn_imputed,
    delimiter=",",
    fmt="%.4f"
)


# Perform PCA
variant_std = gn_imputed.std(axis=1)

print("Zero-variance variants:", np.sum(variant_std == 0))
print("Min std:", variant_std.min())

nonzero_var = gn_imputed.std(axis=1) > 0
gn_pca = gn_imputed[nonzero_var, :]

print(f"Removed {(~nonzero_var).sum()} zero-variance variants")

coordinates, model = allel.pca(gn_pca, n_components=10, scaler='patterson')


# Store PCA-transformed data in a DataFrame (equivalent to pca variable in the original code)
pca = pd.DataFrame(coordinates, columns=[f'PC{i+1}' for i in range(10)])

pca['sample_id'] = samples_filtered
pca.to_csv("pca_unprocessed_projects.csv",index=False)

# Function to clean sample names by removing unwanted patterns
def clean_sample_name(sample_name):
    # Use regex to remove patterns like ../../per_samples/ and '_SXX'
    cleaned_name = re.sub(r".*\/|\'|_S\d+", "", sample_name)
    return cleaned_name


# Access explained variance ratio
explained_variance_ratio = model.explained_variance_ratio_

# Calculate the cumulative variance explained by the components
cumulative_explained_variance = explained_variance_ratio.cumsum()

# Print how much percentage is explained by the first 3 components
print("Explained variance by each component:", explained_variance_ratio)
print("Cumulative variance explained by components:", cumulative_explained_variance)

# Sample names should be included with coordinates, assuming `samples` contains the sample names
sample_names = callset['samples'] # assuming sample names are in g object
gene = callset['variants/GENE']
# Include sample names in coordinates

coordinates_df = pd.DataFrame(coordinates, columns=[f'PC{i+1}' for i in range(10)])
coordinates_df['sample_id'] = samples_filtered
# Apply this cleaning function to your sample names
coordinates_df['sample_id'] = coordinates_df['sample_id'].apply(clean_sample_name)
print(coordinates_df.head())
# Define a function to detect outliers using the IQR method
def detect_outliers_iqr(data):
    Q1 = np.percentile(data, 25, axis=0)
    Q3 = np.percentile(data, 75, axis=0)
    IQR = Q3 - Q1
    lower_bound = Q1 - 2 * IQR
    upper_bound = Q3 + 2 * IQR
    return (data < lower_bound) | (data > upper_bound)

# Detect outliers across all 4 components
outliers_mask = detect_outliers_iqr(coordinates[:, :4])  # Checking the first 4 PCs

# Remove outliers
non_outliers_mask = ~outliers_mask.any(axis=1)
# Sample name to be removed
#sample_to_remove = "ACH-21-FB-067"

# Filter out the sample
#coordinates_no_outliers = coordinates_df[coordinates_df['sample_id'] != sample_to_remove]

coordinates_no_outliers = coordinates_df[non_outliers_mask]

#Print outliers
outlier_indices = np.where(outliers_mask.any(axis=1))[0]
outliers_df = coordinates_df.iloc[outlier_indices]
print("Outliers filtered out:")
print(outliers_df)


# Categorize the project of the samples based on the start of the sample name
def categorize_project(sample_id):
    if sample_id.startswith('ACH21') | sample_id.startswith('ACH-21') | sample_id.startswith('GMTP') | sample_id.startswith('TG') | sample_id.startswith('ACH22'):
        return 'EMAGEN'
    elif sample_id.startswith('ERR1106'):
        return 'Ethiopia_2015_MalariaGEN'
    elif sample_id.startswith('ERR1818'):
        return 'Kenya_2014_MalariaGEN'
    elif sample_id.startswith('ERR2299'):
        return 'Sudan_2015_MalariaGEN'
    elif sample_id.startswith('ERR289'):
        return 'Tanzania_2014_MalariaGEN'
    elif sample_id.startswith('ERR34863'):
        return 'Kenya_2018_MalariaGEN'
    elif sample_id.startswith('CW_'):
        return 'Kenya_2018'
    elif sample_id.startswith('ERR701'):
        return 'Kenya_2014_MalariaGEN'
    elif sample_id.startswith('KU14') | sample_id.startswith('Pf0') | sample_id.startswith('K69'):
        return 'Kenya_2015'
    else:
        return 'Other'

coordinates_no_outliers['site'] = coordinates_no_outliers['sample_id'].apply(categorize_project)


# Filter for samples categorized as 'Other'
samples_in_other = coordinates_no_outliers[coordinates_no_outliers['site'] == 'Other']

# Print or save the sample names
print("Samples in 'Other':")
print(samples_in_other['sample_id'])


coordinates_no_outliers['collectfion_year'] = coordinates_no_outliers['sample_id'].apply(lambda x: 2022 if 'ACH22' in x or 'GMTP' in x else 2021)
#distinct_colors = ['#0000FF', '#FF0000', '#00FF00',  # blue, red, green
                   #'#FFFF00', '#800000', '#808000','#00FFFF',  # yellow, maroon, olive, cyan
                   ##'#FF00FF', '#E7BA52', '#800080', '#808080',  # magenta, teal, purple, gray
                   #'#339966', '#DE9Ed6', '#6B6ECF']             # orange, brown, blueviolet


distinct_colors = {
    'EMAGEN': '#0000FF',                    # blue
    'Kenya_2018': '#FF0000',  # red
    'Ethiopia_2015_MalariaGEN': '#00FF00',     # green
    'Kenya_2014_MalariaGEN': '#FFFF00',     # yellow
    'Sudan_2015_MalariaGEN': '#800000',  # maroon
    'Tanzania_2014_MalariaGEN': '#808000',     # olive
    'Kenya_2018_MalariaGEN': '#00FFFF',                # cyan
    'Kenya_2015': '#FF00FF',                # magenta
    'Other': '#808080', 
}

# Print number of missing positions per site
n_variants = g.shape[0]  # Total number of variants

sample_missingness = pd.DataFrame({
    'sample_id': [clean_sample_name(s) for s in callset['samples'][sample_mask]],
    'missing_count': (missing_per_sample[sample_mask] * n_variants).astype(int),
    'missing_rate': missing_per_sample[sample_mask]
})

sample_missingness['site'] = sample_missingness['sample_id'].apply(categorize_project)

missingness_by_site = sample_missingness.groupby('site').agg({
    'missing_count': 'mean',
    'missing_rate': 'mean',
    'sample_id': 'count'
}).rename(columns={'sample_id': 'n_samples'})

print(f"\n=== Mean Missingness Per Site (out of {n_variants} variants) ===")
print(missingness_by_site.sort_values('missing_count', ascending=False))





# Set up a plot for the first 2 principal components without outliers
# Create the figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(9, 4.0), sharey=False)

# Define colors for each site
colors = {}
unique_sites = coordinates_no_outliers['site'].unique()
#for i, site in enumerate(unique_sites):
    #colors[site] = distinct_colors[i % len(distinct_colors)]

# Plot PC1 vs PC2  in the first subplot
for site in unique_sites:
    subset = coordinates_no_outliers[coordinates_no_outliers['site'] == site]
    axes[0].scatter(subset['PC1'], subset['PC2'], s=15, label=site, alpha=0.5, color=distinct_colors.get(site, '#808080'))

axes[0].set_xlabel(f'PC1 ({explained_variance_ratio[0]*100:.2f}% )')
axes[0].set_ylabel(f'PC2 ({explained_variance_ratio[1]*100:.2f}% )')
#axes[0].set_title("A")

# Plot PC3 vs PC2 in the second subplot
for site in unique_sites:
    subset = coordinates_no_outliers[coordinates_no_outliers['site'] == site]
    axes[1].scatter(subset['PC2'], subset['PC3'], s=15, label=site, alpha=0.5, color=distinct_colors.get(site, '#808080'))

axes[1].set_xlabel(f'PC2 ({explained_variance_ratio[1]*100:.2f}% )')
axes[1].set_ylabel(f'PC3 ({explained_variance_ratio[2]*100:.2f}% )')
#axes[1].set_title("B")

# Place a single legend outside the plots
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, title='Site', markerscale=2,bbox_to_anchor=(0.879, 0.55), loc='center left')

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('PCA_project.png', dpi=300, bbox_inches='tight')
plt.show()


# Print all site-to-color mappings
for site, color in colors.items():
    print(f"Site: {site}, Color: {color}")

