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
        tic;AgNet = panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out,...
            save_temp, alpha, save_pairs, modeProcess);toc;
        
        % Load the expected result
        ExpAgNet = textread('test_data/panda.test.txt');
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        ExpAgNet = reshape(ExpAgNet,[size(AgNet,2), size(AgNet,1)])';% Only in octave

        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(AgNet-ExpAgNet)));
	    assertTrue( deltaMat < tolMat );
        
        % Compare distances in single and double precision
        alpha = 0.7;
        distances={'euclidean','seuclidean','squaredeuclidean','Tfunction','cityblock','minkowski',...
        'chebychev','cosine','correlation',@TfunctionDist};
        % not testing for jaccard, hamming, spearman because for discrete
        % variables
        figure;j=0;
        for distance = distances
            j=j+1;
            distance=distance{1};
            tic;AgNet2 = panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out,...
                    save_temp, alpha, save_pairs, modeProcess,0.5,0,distance,'cpu','single');toc;
            tic;AgNet3 = panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out,...
                    save_temp, alpha, save_pairs, modeProcess,0.5,0,distance,'cpu','double');toc;
            tolMat=3e-5;
            deltaMat=max(max(abs(AgNet2-AgNet3)));
            assert( deltaMat < tolMat );
        end
        
        % Check that gpu implementation gives the same results as CPU
        % implementation
        alpha = 0.7;
        distances={'euclidean','seuclidean','squaredeuclidean','Tfunction','cityblock','minkowski',...
        'chebychev','cosine','correlation','hamming','jaccard','spearman',@TfunctionDist};
        for distance = distances
            distance=distance{1};
            tic;AgNet2 = panda_run(lib_path,exp_file, motif_file, ppi_file, panda_out,...
                    save_temp, alpha, save_pairs, modeProcess,0.5,0,distance,'cpu','double');toc;
            [Exp,RegNet,TFCoop,~,~]=processData(exp_file,motif_file,ppi_file,'union');
            GeneCoReg = Coexpression(Exp); 
            RegNet    = NormalizeNetwork(RegNet);
            GeneCoReg = NormalizeNetwork(GeneCoReg);
            TFCoop  = NormalizeNetwork(TFCoop);
            verbose = 0;
            AgNet3  = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, 0.5, distance, 'cpu', 'double', verbose);
            tolMat  = 1e-14;
            deltaMat= max(max(abs(AgNet2-AgNet3)))
            assert( deltaMat < tolMat );
        end
        
end
