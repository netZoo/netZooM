function test_suite=testPairs2Mat()
        initTestSuite;
end

function testPairs2MatSimple()
        % Add path
        addpath(genpath(fullfile(pwd,'tests')));
        
        pkg load nan;
        nGenes=1000;
        
        % load test panda network
        networkPair='panda.test.pairs.txt_FinalNetwork.pairs';
        priorNetTest=Pairs2Mat(networkPair,nGenes,0);
        
        % Compare networks
        load('panda.test.mat');
        tolMat=1e-6;
        deltaMat=max(max(abs(priorNet-priorNetTest)));
	    assertTrue( deltaMat < tolMat );
end
