function test_suite=testProcessData()
	initTestSuite;
end

function testProcessDataSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
            %we need the nan package because it has a fast implementation of corrcoeff
            %pkg load statistics
            pkg load nan;
            return % readtable is not available in octave
        end

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));
        
        % read expression data
        exp_file   = 'test_data/pantests/expression_down.txt';
        motif_file = 'test_data/pantests/motif.txt';
        ppi_file   = 'test_data/pantests/ppi.txt';
        
        % Intersection
        modeProcess= 'intersection';
        [Exp,~,~,~,GeneNames]=processData(exp_file,motif_file,ppi_file,modeProcess);
        assert(Exp(2,4)==1)
        assert(isequal(GeneNames{2},'gene2'))
        
        % Union
        modeProcess= 'union';
        [Exp,~,~,~,GeneNames]=processData(exp_file,motif_file,ppi_file,modeProcess);
        assert(Exp(2,4)==1)
        assert(isequal(GeneNames{2},'gene2'))
     
end
