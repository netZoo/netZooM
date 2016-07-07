#PUMA
##PANDA Using MicroRNA Associations
PUMA, or **P**ANDA **U**sing **M**icroRNA **A**ssociations, is an extension of the gene regulatory network reconstruction algorithm PANDA, which was published in "Passing messages between biological networks to refine predicted interactions" by Glass K, Huttenhower C, Quackenbush J, Yuan GC. *PLoS One*. 2013 May 31l8(5):e64983, doi: 10.1371/journal.pone.0064832, PMID: 23741402

PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators of mRNA expression levels. This Github repository contains both a C++ version and a MATLAB version of PUMA. The C++ version is based on PANDA version 2. The original PANDA C++ code is available at: http://sourceforge.net/projects/panda-net/. The MATLAB script to run PUMA has fewer options than the C++ code, but runs faster.

##C++ code
The C++ code of PUMA and example files to run PUMA can be found in folder `PUMAc`. The code can be compiled using:
```
g++ PUMA.c -O3 -o PUMA
```

To run PUMA with the assumption that all regulators can form complexes (estimate *responsibility* for all regulators, *eg* a gene regulatory prior with transcription factors only):
```
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -o ToyOutput
```

To tell PUMA to discriminate between regulators that can, and regulators that cannot cannot form complexes (for example a list of microRNAs in `miRlist.txt`), run:
```
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -u miRlist.txt -o ToyOutput
```

###MATLAB code
The MATLAB code of PUMA and example files to run PUMA in MATLAB can be found in folder `PUMAm`. The main script to run PUMA is `RunPUMA.m`. This script is set-up to run PUMA on example data from folder `ToyData`. The script must be run with a regulatory prior (variable `motif_file` in the RunPUMA script) and expression data (variable `exp_file`). Please note that, in principle, PUMA can be run on an identity matrix of expression data as well, but this is not implemented. Protein-protein interaction data (variable `ppi_file`) and a list of microRNAs (variable `mir_file`) are optional parameters. If no `mir_file` is listed, the script will run PANDA to estimate regulatory edges. If a `mir_file` is listed, the script will run PUMA to estimate regulatory edges. In particular, regulators listed in that file will be treated as regulators that cannot form complexes with other regulators, such as microRNAs, while regulators not listed will be treated as regulators that can form complexes, such as transcription factors.

To run PUMA on your own data, change the paths and filenames under "Set Program Parameters". MATLAB functions `NormalizeNetwork.m`, `PANDA.m`, `PUMA.m`, `Tfunction.m`, and `UpdateDiagonal.m` will be called from the main script.

Some important notes of this script compared to the PUMA C++ code:
- The target genes in the regulatory prior (`motif_file`) should match the target genes in the expression data (`exp_file`).
- The regulators in `mir_list` should match symbols in the regulatory prior (`motif_file`).
- Self-interactions in the protein-protein interaction data will be set to 1. The algorithm converges around these edges. Protein-protein interaction "prior" edges therefore should have weights between 0-1.
- In the protein-protein interaction prior, edges between microRNAs from the `mir_file` (see `RunPUMA.m` script) and other regulators can have values, but these will automatically be set to 0, as the algorithm assumes that any regulator listed in `mir_file` will not be able to form edges with other regulators.
