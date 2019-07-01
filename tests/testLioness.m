function test_suite=testLioness()
	initTestSuite;
end

function testLionessSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
		pkg load nan;
        end

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));
        
        % Set Program Parameters
        % exp, motif, and ppi are produced by panda_run with save_temp
        % parameter set != ''
        exp_file   = 'test_data/expression.transposed.mat';
        motif_file = 'test_data/motif.normalized.mat';
        ppi_file   = 'test_data/ppi.normalized.mat';
        panda_file = 'panda2.test.mat'; % This test network has been transposed
        % to test savePairs, so we need to tranpose it back
        load('test_data/panda.test.mat');
        AgNet      = AgNet';
        save('panda2.test.mat','AgNet');
        alpha      = 0.1;
        START      = 1;  % sample-of-interest starting from this index
        END        = 1;  % sample-of-interest ending to this index; use -1 to end at the last sample
        ascii_out  = 0;  % set to 1 if you prefer text output file instead of MAT-file
        save_dir   = 'test_data';
        lib_path   = '../netZooM';

        % Create save folder
        %mkdir tmp;

        % Call Lioness
        lioness_run(exp_file, motif_file, ppi_file, panda_file, save_dir, START, END, alpha, ascii_out, lib_path);

        % Load the computed results
        result = load('lioness.1.mat');
        
        % Load the expected result
        load('test_data/lioness.test.mat');

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(PredNet-result.PredNet)));
	    assertTrue( deltaMat < tolMat );
end
