function test_suite=testPanda()
	initTestSuite;
end

function testPandaSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
		pkg load nan;
        end

        % Set Program Parameters
        exp_file   = 'test_data/expression.txt';
        motif_file = 'test_data/motifTest.txt';
        ppi_file   = 'test_data/ppi.txt';
        panda_out  = '';  % optional, leave empty if file output is not required
        save_temp  = '';  % optional, leave empty if temp data files will not be needed afterward
        lib_path   = '../netZooM';  % path to the folder of PANDA source code
        alpha      = 0.1;
        save_pairs = 0;%saving in .pairs format
        modeProcess= 'intersection';

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));

        % Create save folder
        %mkdir tmp;

        % Call Panda
        AgNet = panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out,...
            save_temp, alpha, save_pairs, modeProcess);

        % Load the expected result
        ExpAgNet = textread('test_data/panda.test.txt');
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet-ExpAgNet)));
	    assertTrue( deltaMat < tolMat );
end
