% Set Program Parameters
exp_file   = '../LIONESS_ENS_GTEX_661_30266/expression.txt';
motif_file = '../LIONESS_ENS_GTEX_661_30266/motif.txt';
ppi_file   = '../LIONESS_ENS_GTEX_661_30266/ppi2015_freeze.txt';
panda_out  = '/tmp/panda.GTEx.mat';  % optional, leave empty if file output is not required
save_temp  = '/tmp';  % optional, leave empty if temp data files are not needed afterward
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;
