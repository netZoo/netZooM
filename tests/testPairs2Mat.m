function test_suite=testPairs2Mat()
        initTestSuite;
end

function testPairs2MatSimple()
        pkg load nan;
        nGenes=1000;
        % load cell line panda network of lymphoblastoid cell lines (LCLs)
        % Reference: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4111-x
        networkPair1='test_data/panda.test.pairs.txt_FinalNetwork.pairs';
        priorNetTest=Pairs2Mat(networkPair1,nGenes,0);
        
        % Compare networks
        load('panda.test.mat');
        tolMat=1e-6;
        deltaMat=max(max(abs(priorNet-priorNetTest)));
	    assertTrue( deltaMat < tolMat );
end
