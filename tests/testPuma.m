function test_suite=testPuma()
	initTestSuite;
end

function testPumaSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
            %we need the nan package because it has a fast implementation of corrcoeff
            %pkg load statistics
            pkg load nan;
        end

        % Set Program Parameters
        outtag='PUMA';
        % Set Program Parameters
        alpha=0.1;
        motif_file='PUMA_ToyMotifData.txt';
        exp_file='PUMA_ToyExpressionData.txt';
        ppi_file='PUMA_ToyPPI.txt';
        % ppi_file=''; % PUMA/PANDA can be run without PPI data
        % mir_file=''; % this will run PANDA
        mir_file='PUMA_ToyMiRList.txt'; % run PUMA

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));

        % Create save folder
        %mkdir tmp;

        % Call Puma
        AgNet=RunPUMA(outtag,alpha,motif_file,exp_file,ppi_file,mir_file);

        % Load the expected result
        ExpAgNet = textread('puma.test.txt');
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet-ExpAgNet)));
	    assertTrue( deltaMat < tolMat );
end

function testPumaLioness()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
            pkg load nan;
        end

        % Set Program Parameters
        outtag='Test';
        % Set Program Parameters
        alpha=0.1;
        motif_file='PUMA_ToyMotifData.txt';
        exp_file='PUMA_ToyExpressionData.txt';
        ppi_file='PUMA_ToyPPI.txt';
        % ppi_file=''; % PUMA/PANDA can be run without PPI data
        % mir_file=''; % this will run PANDA
        mir_file='PUMA_ToyMiRList.txt'; % run PUMA

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));

        % Create save folder
        %mkdir tmp;

        % Call Puma lioness
        RunPUMALIONESS(outtag,alpha,motif_file,exp_file,ppi_file,mir_file);

        % Load the expected result
        AgNet    = load('PUMA_LIONESSNetworks.mat');
        ExpAgNet = load('Test_LIONESSNetworks.mat');
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        %ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet.PredNet-ExpAgNet.PredNet)));
	    assertTrue( deltaMat < tolMat );
end

function testPumaLionessSubset()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
            pkg load nan;
        end

        % Set Program Parameters
        outtag='PUMA2';
        % Set Program Parameters
        alpha=0.1;
        motif_file='PUMA_ToyMotifData.txt';
        exp_file='PUMA_ToyExpressionData.txt';
        ppi_file='PUMA_ToyPPI.txt';
        % ppi_file=''; % PUMA/PANDA can be run without PPI data
        % mir_file=''; % this will run PANDA
        mir_file='PUMA_ToyMiRList.txt'; % run PUMA
        SelectSize=20; % nr of samples to run
        Offset=0; % offset

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));

        % Create save folder
        %mkdir tmp;

        % Call Puma lioness
        RunPUMALIONESSsubset(outtag, alpha, motif_file, exp_file, ppi_file, mir_file, SelectSize, Offset);

        % Load the expected result
        AgNet    = load('PUMA2_LIONESSNetworks.mat');
        ExpAgNet = load('Test2_LIONESSNetworks.mat');
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        %ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet.PredNet-ExpAgNet.PredNet)));
	    assertTrue( deltaMat < tolMat );
end

% Commented testPumaLionessRs because octave does not have rng
%function testPumaLionessRs()
%	% Tell if this is Octave (Unit tests) or Matlab
%        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%        % Load statistics package from Octave
%        if isOctave
%                %we need the nan package because it has a fast implementation of corrcoeff
%                %pkg load statistics
%            pkg load nan;
%        end

%        % Set Program Parameters
%        outtag='PUMA_rs';
%        % Set Program Parameters
%        alpha=0.1;
%        motif_file='PUMA_ToyMotifData.txt';
%        exp_file='PUMA_ToyExpressionData.txt';
%        ppi_file='PUMA_ToyPPI.txt';
%        % ppi_file=''; % PUMA/PANDA can be run without PPI data
%        % mir_file=''; % this will run PANDA
%        mir_file='PUMA_ToyMiRList.txt'; % run PUMA
%        perc=10; % percentage of samples to remove
%        nrboots=10; % how many times to resample the data

%        % Add path
%        addpath(genpath(fullfile(pwd,'tests')));

%        % Create save folder
%        %mkdir tmp;

%        % Call Puma lioness
%        rng(420)%set the seed
%        PredNet=RunPUMAresample(outtag, alpha, motif_file, exp_file, ppi_file, mir_file, perc, nrboots);

%        % Load the expected result
%        ExpAgNet = load('PUMA_rs_test.mat');
%        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
%        %ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';

%        % Compare the outputs
%        tolMat=1e-6;
%        deltaMat=max(max(abs(PredNet-ExpAgNet.PredNet)));
%	    assertTrue( deltaMat < tolMat );
%end