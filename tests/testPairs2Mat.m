function test_suite=testPairs2Mat()
        initTestSuite;
end

function testPairs2MatSimple()
        pkg load nan;
        nGenes=1000;
        % load test panda network
        networkPair1='test_data/panda.test.pairs.txt_FinalNetwork.pairs';
        priorNetTest=Pairs2Mat(networkPair1,nGenes,0);
        
        % Compare networks
        load('panda.test.mat');
        tolMat=1e-6;
        deltaMat=max(max(abs(priorNet-priorNetTest)));
	    assertTrue( deltaMat < tolMat );
end
