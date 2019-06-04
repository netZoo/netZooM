% Tell if this is Octave (Unit tests) or Matlab
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Load statistics package from Octave
if isOctave
	pkg load statistics
end

% Set Program Parameters
exp_file   = 'test_data/expression.txt';
motif_file = 'test_data/motif.txt';
ppi_file   = 'test_data/ppi.txt';
panda_out  = 'tmp/panda.test.txt';  % optional, leave empty if file output is not required
save_temp  = '';  % optional, leave empty if temp data files will not be needed afterward
lib_path   = '../netzoo-m';  % path to the folder of PANDA source code
alpha      = 0.1;
save_pairs = 0;%saving in .pairs format

% Add path
addpath(genpath(fullfile(pwd,'../netzoo-m')));

% Create save folder
pwd
mkdir tmp;

% Call Panda
AgNet=panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out, save_temp, alpha, save_pairs);
