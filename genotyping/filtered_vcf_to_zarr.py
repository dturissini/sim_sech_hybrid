# Import packages.
import allel
import numcodecs
import numpy as np
import sys
import zarr


# Define file paths.
chrom = str(sys.argv[1])
vcf_path = str(sys.argv[2])
zarr_path = str(sys.argv[3])


# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(
    vcf_path, zarr_path, group=str(chrom),
    fields=['samples', 'DP', 'POS', 'CHROM', 'GT', 'variants/REF', 'variants/ALT'], log=sys.stdout, overwrite=True,
)
