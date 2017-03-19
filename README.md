# MATLAB PANDA

MATLAB implementation of PANDA & LIONESS algorithms.

Full PANDA & LIONESS algorithms were described in the following literature:

* Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PloS one 8.5 (2013): e64832.
* Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." arXiv preprint arXiv:1505.06440 (2015).

Author: cychen (ntu.joey@gmail.com), marieke, kimbie.

Orignial source code adapted from marieke & kimbie's version.

# Usage

## PANDA

1. Set up PANDA run-time parameters by editing `panda_config.m`.
2. Run PANDA main program via `panda_run.sh`.

## LIONESS

1. Run PANDA first to get preprocessed middle files and aggregated PANDA network.
2. Set up LIONESS run-time parameters by editing `lioness_config.m`.
3. Run LIONESS main program via `lioness_run.sh`.

# File format

See example input files in `test_data/`.

The example files include following data:

1. A toy expression profile (1000 genes x 50 samples).
2. A list of genome-wide TF-target interactions (motif prior). 
3. A list of protein-protein interactions (PPIs) between TFs.
4. The output PANDA network building on the example data.

# What's new in this version

We reach a 20-25% reduction in computation time and memory usage in comparison to its early version. The following optimazation has been implemented in this version:

## Network Normalization

Use MATLAB built-in zscore function instead of calculation from scratch. -> 30% speed-up.

## T-function

Use bsxfun instead of repmat, using symmetric matrix multiplication, and reusing summed-square vector. -> 25% speed-up.

## PANDA Function

Move out network normalization from PANDA function to main program to remove unnecessary repeated normalization in the following LIONESS-loop: both PPI network and motif network need to be normalized only once. -> 1-10 sec reduced for each LIONESS iteration depending on the network sizes.

Save W matrix (R+A) so that we don't have to compute it twice. -> ~0.5% overall speed-up.

Check hamming inside the while-loop to reduce the last unnecessary updates on TFcoop and GeneCoReg networks. -> ~2% overall speed-up.

## I/O

Use MATLAB binary files (mat-file) instead of text files for I/O can boost the performance (>10x faster).

## Reusing the processed middle files

Save the input expression matrix (transposed), normalized motif/PPI networks, and the aggregated PANDA network to binary files (mat-file) for later use in each LIONESS run: expression matrix can only be saved to v7.3 mat file (compressed Matlab HDF5 format) when its size is over 2GB; normalized motif/PPI networks and output PANDA network can be saved as v6 uncompressed mat files for faster I/O efficiency.

## Computing gene-gene coexpression matrix

Move out this part as an independent function so that no need to transpose the input matrix each time in LIONESS iteration. Consequently, we save a copy of expression matrix in the memory and also one matrix transpose operation overhead for each single-sample LIONESS run.
