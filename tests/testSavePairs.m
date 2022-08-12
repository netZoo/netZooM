function test_suite=testSavePairs()
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function testSavePairsSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
            pkg load nan;
        end
        %1st case
        TFNames  ={'a','b'};
        GeneNames={'l','m'};
        AgNet    =magic(2);
        RegNet   =[];
        outtag   ='';
        [TF,gene,AgNetCol,RegNetCol]=SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag);
        deltaMat=max(max(AgNetCol-AgNet(:)));
        assert(deltaMat<1e-6);
        assert(all(strcmpi(TF,{'a','b','a','b'}')));
        assert(all(strcmpi(gene,{'l','l','m','m'}')));
        assert(isempty(RegNetCol));
        %2nd case
        TFNames  ={'a','b'};
        GeneNames={'l','m'}';
        [TF,gene,AgNetCol,RegNetCol]=SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag);
                deltaMat=max(max(AgNetCol-AgNet(:)));
        assert(deltaMat<1e-6);
        assert(all(strcmpi(TF,{'a','b','a','b'}')));
        assert(all(strcmpi(gene,{'l','l','m','m'}')));
        assert(isempty(RegNetCol));
        %3rd case
        TFNames  ={'a','b'}';
        GeneNames={'l','m'}';
        [TF,gene,AgNetCol,RegNetCol]=SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag);
                deltaMat=max(max(AgNetCol-AgNet(:)));
        assert(deltaMat<1e-6);
        assert(all(strcmpi(TF,{'a','b','a','b'}')));
        assert(all(strcmpi(gene,{'l','l','m','m'}')));
        assert(isempty(RegNetCol));
        %4th case
        TFNames  ={'a','b'}';
        GeneNames={'l','m'};
        [TF,gene,AgNetCol,RegNetCol]=SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag);
                deltaMat=max(max(AgNetCol-AgNet(:)));
        assert(deltaMat<1e-6);
        assert(all(strcmpi(TF,{'a','b','a','b'}')));
        assert(all(strcmpi(gene,{'l','l','m','m'}')));
        assert(isempty(RegNetCol));
end
