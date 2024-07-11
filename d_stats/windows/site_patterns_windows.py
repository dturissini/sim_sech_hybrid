# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[0] = reference
# sys.argv[1] = p1 population
# sys.argv[2] = p2 population
# sys.argv[3] = p3 population
# sys.argv[4] = p4 population
# sys.argv[5] = zarr_file
# sys.argv[6] = windows_file
# sys.argv[7] = metadata_file

# Define a function to load genotyope and positions arrays.
def load_callset_pos(chrom, zarr_file):
    # Load the vcf file.
    callset = zarr.open_group(zarr_file, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos


# Define a function to calculate alternative allele frequencies.
def calc_alt_freqs(gt):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs

# Define a site pattern function.
def site_patterns(p1, p2, p3, p4):
    # Calculate site pattern counts.
    abba = np.nansum((1 - p1) * (p2) * (p3) * (1 - p4))
    baba = np.nansum((p1) * (1 - p2) * (p3) * (1 - p4))
    baaa = np.nansum((p1) * (1 - p2) * (1 - p3) * (1 - p4))
    abaa = np.nansum((1 - p1) * (p2) * (1 - p3) * (1 - p4))
    return abba, baba, baaa, abaa

# Define a function to calculate site patterns.
def dros_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx, p4_idx,
):
    # Determine the indicies where each population has called genotypes.
    p1_mask = (gt.take(p1_idx, axis=1).is_called() == True).any(axis=1)
    p2_mask = (gt.take(p2_idx, axis=1).is_called() == True).any(axis=1)
    p3_mask = (gt.take(p3_idx, axis=1).is_called() == True).any(axis=1)
    p4_mask = (gt.take(p4_idx, axis=1).is_called() == True).any(axis=1)
    # Determine the indicied where all populations have called genotypes.
    called_mask = (p1_mask & p2_mask & p3_mask & p4_mask)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to np.nan since we don't have any sites to perform computations on.
        results = np.full(4, np.nan)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(4)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p4_alt_freqs = calc_alt_freqs(gt.take(p4_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the allele frequencies based on the most common allele in the outgroup.
            p1_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            p4_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p4_alt_freqs - 1), p4_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, baaa, abaa = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs, p4_der_freqs,
            )
            # Intialize a results array.
            results = np.array([abba, baba, baaa, abaa])
    return results


# Extract the reference assembly and focal populations.
p1 = str(sys.argv[1])
p2 = str(sys.argv[2])
p3 = str(sys.argv[3])
p4 = str(sys.argv[4])
zarr_prefix = str(sys.argv[5])
win_file = str(sys.argv[6])
meta_file = str(sys.argv[7])
results_file = str(sys.argv[8])

# Read in meta data has a pandas dataframe.
meta_df = pd.read_csv(meta_file, sep='\t')
# Intialize a dictionary.
idx_dicc = {}
# For all populations.
for pop in [p1, p2, p3, p4]:
    # Fill the dictiasdfonary.
    idx_dicc[pop] = meta_df[meta_df['population'] == pop].index.values
### Fill the dictionary for the last population.
##idx_dicc['yak'] = np.concatenate((idx_dicc['yak_allo'], idx_dicc['yak_symp']))


# Read in the ortholog dataframe.
window_df = pd.read_csv(win_file, sep='\t')
# Extract ortholog information.
chroms = window_df.chr.values
starts = window_df.start.values
ends = window_df.end.values
# Intialize a results matrix to store the results.
results_mat = np.empty((window_df.shape[0], 6))

# For every ortholog window.
for idx in range(window_df.shape[0]):
    # Extract the ortholog information.
    chrom = chroms[idx]
    start = starts[idx]
    end = ends[idx]
    # Extract the genotype callset and positions.
    zarr_file = zarr_prefix + '_' + chrom + '.zarr'
    callset, all_pos = load_callset_pos(chrom, zarr_file)
    # Identify the window to extract.
    wind_loc = all_pos.locate_range(start, end)
    # Calculate site patterns.
    sps = dros_site_patterns(
        gt=allel.GenotypeArray(callset[wind_loc]),
        p1_idx=idx_dicc[p1], p2_idx=idx_dicc[p2],
        p3_idx=idx_dicc[p3], p4_idx=idx_dicc[p4],
    )
    # Append the results.
    results_mat[idx, :] = sps
        
# Export the the results matrix.
np.savetxt(
    results_file,
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)