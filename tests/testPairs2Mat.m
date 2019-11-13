function test_suite=testPairs2Mat()
        initTestSuite;
end

function testPairs2MatSimple()
        % Add path
        addpath(genpath(fullfile(pwd,'tests')));
        
        % Load packages
        if isOctave
            pkg load nan;
        end
        nGenes=1000;
        
        % load test panda network
        networkPair='panda.test.pairs.txt_FinalNetwork.pairs';
        priorNetTest=Pairs2Mat(networkPair,nGenes,0);
        
        % Compare networks
        load('panda.test.mat');
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet-priorNetTest)));
	    assertTrue( deltaMat < tolMat );
end
