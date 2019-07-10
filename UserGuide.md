# MATLAB PANDA

MATLAB implementation of PANDA & LIONESS algorithms.

Full PANDA & LIONESS algorithms were described in the following literature:

* Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PLoS ONE 8.5 (2013): e64832.
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

Use bsxfun instead of repmat, use symmetric matrix multiplication, and reuse summed-square vector. -> 25% speed-up.

## PANDA Function

Move out network normalization from PANDA function to main program to remove unnecessary repeated normalization in each following LIONESS iteration: both PPI network and motif network need to be normalized only once. -> 1-10 sec reduced for each LIONESS iteration depending on the network sizes.

Save W matrix (R+A) so that we don't have to compute it twice. -> ~0.5% overall speed-up.

Check hamming inside the while-loop to reduce the last unnecessary updates on TFcoop and GeneCoReg networks. -> ~2% overall speed-up.

## I/O

Use MATLAB binary files (mat-file) instead of text files for I/O which boost the performance (>10x faster).

## Reusing the processed middle files

Save the input expression matrix (transposed), normalized motif/PPI networks, and the aggregated PANDA network to binary files (mat-file) for later use in each LIONESS run: expression matrix can only be saved to v7.3 mat file (compressed Matlab HDF5 format) when its size is over 2GB; normalized motif/PPI networks and output PANDA network can be saved as v6 uncompressed mat files for best I/O efficiency.

## Computing gene-gene coexpression matrix

Move out this part as an independent function so that no need to transpose the input matrix each time in LIONESS iteration. Consequently, we save a copy of expression matrix in the memory and also one matrix transpose operation overhead for each single-sample LIONESS run.


# PUMA #
## PANDA Using MicroRNA Associations ##
PUMA, or **P**ANDA **U**sing **M**icroRNA **A**ssociations, is an extension of the gene regulatory network reconstruction algorithm PANDA, which was published in "Passing messages between biological networks to refine predicted interactions" by Glass K, Huttenhower C, Quackenbush J, Yuan GC. *PLoS One*. 2013 May 31l8(5):e64983, doi: 10.1371/journal.pone.0064832, PMID: 23741402

PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators of mRNA expression levels. This Github repository contains both a C++ version and a MATLAB version of PUMA. The C++ version is based on PANDA version 2. The original PANDA C++ code is available at: http://sourceforge.net/projects/panda-net/. The MATLAB script to run PUMA has fewer options than the C++ code, but runs faster.

## MATLAB code ##
### Condition-specific PUMA networks ###
The MATLAB code of PUMA and example files to run PUMA in MATLAB can be found in folder `PUMAm`. The main script to run PUMA is `RunPUMA.m`. This script is set-up to run PUMA on example data from folder `ToyData`. The script must be run with a regulatory prior (variable `motif_file` in the RunPUMA script) and expression data (variable `exp_file`). Please note that, in principle, PUMA can be run on an identity matrix of expression data as well, but this is not implemented. Protein-protein interaction data (variable `ppi_file`) and a list of microRNAs (variable `mir_file`) are optional parameters. If no `mir_file` is listed, the script will run PANDA to estimate regulatory edges. If a `mir_file` is listed, the script will run PUMA to estimate regulatory edges. In particular, regulators listed in that file will be treated as regulators that cannot form complexes with other regulators, such as microRNAs, while regulators not listed will be treated as regulators that can form complexes, such as transcription factors.

To run PUMA on your own data, change the paths and filenames under "Set Program Parameters". MATLAB functions `NormalizeNetwork.m`, `PANDA.m`, `PUMA.m`, `Tfunction.m`, and `UpdateDiagonal.m` will be called from the main script.

Some important notes of this script compared to the PUMA C++ code:
- The target genes in the regulatory prior (`motif_file`) should match the target genes in the expression data (`exp_file`).
- The regulators in `mir_list` should match symbols in the regulatory prior (`motif_file`).
- Self-interactions in the protein-protein interaction data will be set to 1. The algorithm converges around these edges. Protein-protein interaction "prior" edges therefore should have weights between 0-1.
- In the protein-protein interaction prior, edges between microRNAs from the `mir_file` (see `RunPUMA.m` script) and other regulators can have values, but these will automatically be set to 0, as the algorithm assumes that any regulator listed in `mir_file` will not be able to form edges with other regulators.

### Resampling PUMA networks ###
`RunPUMAresample.m` runs PUMA as described above, but resamples the data multiple times (variable `nrboots` under "Set Program Parameters") by removing a certain percentage (variable `perc` under "Set Program Parameters") of samples to obtain a collection of networks.

### Single-sample PUMA networks ###
LIONESS, or **L**inear **I**nterpolation to **O**btain **N**etwork **E**stimates for **S**ingle **S**amples, can be used to estimate single-sample networks using aggregate networks made with any network reconstruction algorithm (http://arxiv.org/pdf/1505.06440.pdf).

The main script to run PUMA+LIONESS is `RunPUMALIONESS.m`. For instructions to run this script, see instructions under "Condition-specific PUMA networks".

### Other scripts ###
Some useful scripts can be found in the folder `OtherScripts`:
- bash scripts `RunPUMA.sh` and `RunPUMALIONESS.sh` can be used to remotely run `RunPUMA.m` and `RunPUMALIONESS.m`, respectively.
- `getCompleteEdgelist.R` is an R script that can convert an unweighted regulatory prior in a complete edgelist. This can be useful when preparing the regulatory prior and expression data.
