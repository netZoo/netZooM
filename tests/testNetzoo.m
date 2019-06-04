% Set Program Parameters
exp_file   = 'test_data/expression.txt';
motif_file = 'test_data/motif.txt';
ppi_file   = 'test_data/ppi.txt';
panda_out  = '/tmp/panda.test.txt';  % optional, leave empty if file output is not required
save_temp  = '/tmp';  % optional, leave empty if temp data files will not needed afterward
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;

save_pairs=0;
AgNet=panda_run(exp_file, motif_file, ppi_file, panda_out, save_temp, alpha, save_pairs);