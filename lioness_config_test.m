% Set Program Parameters
exp_file   = '/tmp/expression.transposed.mat';  % produced by panda_run
motif_file = '/tmp/motif.normalized.mat';  % produced by panda_run
ppi_file   = '/tmp/ppi.normalized.mat';  % produced by panda_run
panda_file = '/tmp/panda.test.mat';  % produced by panda_run
save_dir   = '/tmp/lioness.test/';  % output dir for LIONESS networks
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;
START      = 1;  % sample-of-interest starting from this index
END        = 5;  % sample-of-interest ending to this index; use -1 to end at the last sample
ascii_out  = 0;  % set to 1 if you prefer text output file instead of MAT-file
