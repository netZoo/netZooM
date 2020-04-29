function test_suite=testOtter()
	initTestSuite;
end

function testOtterSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
            %we need the nan package because it has a fast implementation of corrcoeff
            %pkg load statistics
            pkg load nan;
        end

        % Set Program Parameters
        generatePriors=0;
        if generatePriors==1
            exp_file   = 'test_data/expression.txt';
            motif_file = 'test_data/motifTest.txt';
            ppi_file   = 'test_data/ppi.txt';
            [Exp,RegNet,TFCoop,TFNames,GeneNames]=processData(exp_file,motif_file,ppi_file,'intersection');
            tic; GeneCoReg = Coexpression(Exp); toc;
            csvwrite('c.csv',GeneCoReg);
            csvwrite('p.csv',TFCoop);
            csvwrite('w.csv',RegNet);
        elseif generatePriors==0
            C=csvread('tests/test_data/otter/c.csv');
            P=csvread('tests/test_data/otter/p.csv');
            W=csvread('tests/test_data/otter/w.csv');
        end

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));

        % Call Otter
        tic;W = otter(W,P,C);toc;
        
        % Load the expected result
        filename = 'tests/test_data/otter/test_otter.csv';
        W_test = csvread(filename);

        % Compare the outputs
        tolMat  =1e-6;
        deltaMat=max(max(abs(W-W_test)));
	    assertTrue( deltaMat < tolMat );
        
end
