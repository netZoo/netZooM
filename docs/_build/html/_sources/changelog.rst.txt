==========
Changelog
==========

0.5.1 (2020-06-27)
------------------
- Major fix for OTTER behavior across platforms.

0.5.0 (2020-05-02)
------------------
- new tool: OTTER
- (Cross-platform) Unit test for OTTER

0.4.3 (2020-03-06)
------------------

- gpuPANDA and gpuLIONESS
- gpuPANDA and gpuLIONESS tutorials 

0.4.2 (2020-02-21)
------------------

- optPANDA 
- optPANDA tutorial

0.4.1 (2020-01-18)
------------------

- Cross-plateform haramonization of PANDA behavior.
- Data processing was set to 'Union' by default.

0.4.0 (2019-12-10)
------------------

- SPIDER

0.3.0 (2019-11-17)
------------------

- PUMA
- differential targeting tutorial
- GRAND database tutorial

0.2.0 (2019-11-13)
------------------

- LIONESS
- PANDA tutorial
- LIONESS tutorial

0.1.0 (2019-5-28)
------------------

- Changelog added to the doc

- PANDA's CY Chen implementation that is different from the original implementation on the following points

  We reach a 20-25% reduction in computation time and memory usage in comparison to its early version. The following optimazation has been implemented in this version:

  Network Normalization:

  Use MATLAB built-in zscore function instead of calculation from scratch. -> 30% speed-up.

  T-function:

  Use bsxfun instead of repmat, use symmetric matrix multiplication, and reuse summed-square vector. -> 25% speed-up.

  PANDA Function:

  Move out network normalization from PANDA function to main program to remove unnecessary repeated normalization in each following LIONESS iteration: both PPI network and motif network need to be normalized only once. -> 1-10 sec reduced for each LIONESS iteration depending on the network sizes.

  Save W matrix (R+A) so that we don't have to compute it twice. -> ~0.5% overall speed-up.

  Check hamming inside the while-loop to reduce the last unnecessary updates on TFcoop and GeneCoReg networks. -> ~2% overall speed-up.

  I/O:

  Use MATLAB binary files (mat-file) instead of text files for I/O which boost the performance (>10x faster).

  Reusing the processed middle files:

  Save the input expression matrix (transposed), normalized motif/PPI networks, and the aggregated PANDA network to binary files (mat-file) for later use in each LIONESS run: expression matrix can only be saved to v7.3 mat file (compressed Matlab HDF5 format) when its size is over 2GB; normalized motif/PPI networks and output PANDA network can be saved as v6 uncompressed mat files for best I/O efficiency.

  Computing gene-gene coexpression matrix:

  Move out this part as an independent function so that no need to transpose the input matrix each time in LIONESS iteration. Consequently, we save a copy of expression matrix in the memory and also one matrix transpose operation overhead for each single-sample LIONESS run.


