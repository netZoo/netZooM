% Set Program Parameters
exp_file   = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/expression.transposed.mat';  % produced by panda_run
motif_file = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/motif.normalized.mat';  % produced by panda_run
ppi_file   = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/ppi.normalized.mat';  % produced by panda_run
panda_file = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/panda.mat';  % produced by panda_run
save_dir   = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/networks';  % output dir for LIONESS networks
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;
START      = 1;  % sample-of-interest starting from this index
END        = 5;  % sample-of-interest ending to this index; use -1 to end at the last sample
ascii_out  = 0;  % set to 1 if you prefer text output file instead of MAT-file
