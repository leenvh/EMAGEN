
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

parser = argparse.ArgumentParser(description='Analyse VCF file')
parser.add_argument('vcf', type=str, help='VCF file')

#Example
#python /scripts/analyse_vcf_project.py genotyped_merged.ann.vcf.gz

args = parser.parse_args()

callset = allel.read_vcf(args.vcf)

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


vcf = pysam.VariantFile(args.vcf)
ann = {}
for var in vcf:
    for d in  var.info['ANN']:
        d = d.split('|')
        key = (var.chrom,var.pos,var.alts[0])
        if key in ann and ann[key][2] in ['missense_variant', 'synonymous_variant']:
            continue
        ann[key] = (d[3].replace('gene:',''),d[10].replace('p.',''),d[1])


variant_idx = [
    i for i in range(len(callset['variants/POS'])) 
    if ann[(callset['variants/CHROM'][i], callset['variants/POS'][i], callset['variants/ALT'][i][0])][2] in ['missense_variant', 'synonymous_variant','upstream_gene_variant']
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


def _get_LD(callset: dict) -> np.ndarray:
    g = allel.GenotypeArray(callset['calldata/GT'])
    gn = g.to_n_alt(fill=-1)
    r = allel.rogers_huff_r(gn)
    r2 = squareform(r ** 2)
    return r2

def calculate_LD(callset: dict) -> pd.DataFrame:
    g = allel.GenotypeArray(callset['calldata/GT'])
    gn = g.to_n_alt(fill=-1)
    r = allel.rogers_huff_r(gn)
    r2 = squareform(r ** 2)
    # pickle.dump(r2,open("rsquared.pickle","wb"))

    data = {
        "gene1": [],
        "change1": [],
        "gene2": [],
        "change2": [],
        "LD": []
    }

    for i in range(len(r2)):
        for j in range(len(r2)):
            if r2[i,j] >=0:
                data["gene1"].append(callset['variants/GENE'][i])
                data["change1"].append(callset['variants/AA'][i])
                data["gene2"].append(callset['variants/GENE'][j])
                data["change2"].append(callset['variants/AA'][j])
                data["LD"].append(r2[i,j])

    df = pd.DataFrame(data)

    return df


##### PCA 

#print(callset['variants/AA'])
g = allel.GenotypeArray(callset['calldata/GT'])

# Convert genotype data to allele counts
gn = g.to_n_alt()

# Center the data
gn_centered = gn - np.nanmean(gn, axis=0)

# Perform PCA on the centered data
coordinates, model = allel.pca(gn_centered, n_components=4, scaler=None)

# Store PCA-transformed data in a DataFrame (equivalent to pca variable in the original code)
pca = pd.DataFrame(coordinates, columns=[f'PC{i+1}' for i in range(4)])
pca['sample_id'] = callset['samples']



# Perform PCA with 4 components

pca['sample_id'] = callset['samples']
pca.to_csv("pca_unprocessed.csv",index=False)

# Function to clean sample names by removing unwanted patterns
def clean_sample_name(sample_name):
    # Use regex to remove patterns like ../../per_samples/ and '_SXX'
    cleaned_name = re.sub(r".*\/|\'|_S\d+", "", sample_name)
    return cleaned_name




coordinates, model = allel.pca(g.to_n_alt(), n_components=4, scaler='patterson')

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

coordinates_df = pd.DataFrame(coordinates, columns=[f'PC{i+1}' for i in range(4)])
coordinates_df['sample_id'] = sample_names
# Apply this cleaning function to your sample names
coordinates_df['sample_id'] = coordinates_df['sample_id'].apply(clean_sample_name)
print(coordinates_df.head())
# Define a function to detect outliers using the IQR method
def detect_outliers_iqr(data):
    Q1 = np.percentile(data, 25, axis=0)
    Q3 = np.percentile(data, 75, axis=0)
    IQR = Q3 - Q1
    lower_bound = Q1 - 3 * IQR
    upper_bound = Q3 + 3 * IQR
    return (data < lower_bound) | (data > upper_bound)

# Detect outliers across all 4 components
outliers_mask = detect_outliers_iqr(coordinates[:, :4])  # Checking the first 4 PCs

# Remove outliers
non_outliers_mask = ~outliers_mask.any(axis=1)

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



coordinates_no_outliers['collectfion_year'] = coordinates_no_outliers['sample_id'].apply(lambda x: 2022 if 'ACH22' in x or 'GMTP' in x else 2021)
distinct_colors = ['#0000FF', '#FF0000', '#00FF00',  # blue, red, green
                   '#FFFF00', '#800000', '#808000','#00FFFF',  # yellow, maroon, olive, cyan
                   '#FF00FF', '#E7BA52', '#800080', '#808080',  # magenta, teal, purple, gray
                   '#339966', '#DE9Ed6', '#6B6ECF']             # orange, brown, blueviolet




# Set up a plot for the first 2 principal components without outliers
# Create the figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(9, 4.0), sharey=False)

# Define colors for each site
colors = {}
unique_sites = coordinates_no_outliers['site'].unique()
for i, site in enumerate(unique_sites):
    colors[site] = distinct_colors[i % len(distinct_colors)]

# Plot PC1 vs PC2  in the first subplot
for site in unique_sites:
    subset = coordinates_no_outliers[coordinates_no_outliers['site'] == site]
    axes[0].scatter(subset['PC1'], subset['PC2'], s=15, label=site, alpha=0.5, color=colors[site])

axes[0].set_xlabel(f'PC1 ({explained_variance_ratio[0]*100:.2f}% )')
axes[0].set_ylabel(f'PC2 ({explained_variance_ratio[1]*100:.2f}% )')
#axes[0].set_title("A")

# Plot PC3 vs PC2 in the second subplot
for site in unique_sites:
    subset = coordinates_no_outliers[coordinates_no_outliers['site'] == site]
    axes[1].scatter(subset['PC2'], subset['PC3'], s=15, label=site, alpha=0.5, color=colors[site])

axes[1].set_xlabel(f'PC2 ({explained_variance_ratio[1]*100:.2f}% )')
axes[1].set_ylabel(f'PC3 ({explained_variance_ratio[2]*100:.2f}% )')
#axes[1].set_title("B")

# Place a single legend outside the plots
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, title='Project', bbox_to_anchor=(0.879, 0.55), loc='center left')

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('PCA_project.png', dpi=300, bbox_inches='tight')
plt.show()


# Print all site-to-color mappings
for site, color in colors.items():
    print(f"Site: {site}, Color: {color}")



