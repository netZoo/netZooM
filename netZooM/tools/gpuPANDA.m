function RegNet = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric,...
                computing, precision, verbose, saveMemory, gpuSave, isPUMA,...
                s1, s2, t1, t2)
% Description:
%             GPU-accelerated PANDA, slightly different implmentation that is 
%             optimized for memory.
% Inputs:
%             RegNet    : motif prior of gene-TF regulatory network
%             GeneCoReg : gene-gene co-regulatory network
%             TFCoop    : PPI binding between transcription factors
%             respWeight: real number between 0 and 1. Weight of the responsability matrix (default: 0.5)
%             similarityMetric: string containing the similarity metric that PANDA uses to 
%                         find agreement between networks. Similarity
%                         scores are kept as is, and distance scores were
%                         converted to similarities through s=1/(1+d)
%                         'Tfunction'   : (Default) Modified tanimoto
%                                        similarity as described in doi:10.1371/journal.pone.0064832
%                         @TfunctionDist: same as 'Tfunction' but works
%                                        with pdist2, slower.
%                         'euclidean'  : 1/(1+ euclidean distance). GPU enabled
%                         'squaredeuclidean': 1/(1+ squared euclidean
%                                             distance). GPU enabled
%                         'seuclidean' : 1/(1+ standardized euclidean
%                                        distance) . GPU enabled
%                         'cityblock'  : 1/(1+ cityblock distance). GPU enabled
%                         'minkowski'  : 1/(1+ minkowski distance). 
%                                        with p=3. p=1,p=2, and p=inf are
%                                        covered by cityblock, euclidean,
%                                        and Chebychev respectively. GPU enabled
%                         'chebychev'  : 1/(1+ chebychev distance). GPU enabled
%                         'cosine'     : cosine of the included angle.GPU enabled
%                         'correlation': sample correlation between
%                                        points. GPU enabled.
%                         'hamming'    : 1/(1+ hamming distance). GPU enabled
%                         'jaccard'    : Jaccard coefficient. GPU enabled
%                         'spearman'   : sample Spearman's rank correlation
%             computing : 'cpu'(default)
%                         'gpu' uses GPU to compute distances
%             distance  : computing precision
%                         double: double precision(Default)
%                         single: single precision
%             verbose   : 1 prints iterations (Default)
%                         0 does not print iterations
%             saveMemory: 1 saves memory on device but slower (Default)
%                         0 faster computation but more memory required
%             gpuSave   : GPU flag for saving output network
%                         1 keep the final network in GPU memory
%                         0 (Default) send the final network to CPU and reset GPU memory
%             isPUMA    : boolean, if gpuPANDA is called through PUMA then add
%                          additional steps to keep the miRNA interactions.
%             s1,s2,t1,t2: when isPUMA==1, these are the indices of miR 
%                          interactions in TFCoop, the TF-TF PPI network.
% Outputs:
%             RegNet    : inferred gene-TF regulatory network
% Author(s):
%             Marouen Ben Guebila

    if nargin<5
        respWeight=0.5;
    end
    if nargin<6
        similarityMetric='Tfunction';
    end
    if nargin<7
        computing='cpu';
    end
    if nargin<8
        precision='double';
    end
    if nargin<9
        verbose=1;
    end
    if nargin<10
        saveMemory=1;
    end
    if nargin<11
        gpuSave=0;
    end
    if nargin<12
       isPUMA=0; 
    end
    if iscategorical(similarityMetric)
        similarityMetric=char(similarityMetric(1));
    end
    similarityMetricChar=similarityMetric;
    if isa(similarityMetric,'function_handle')
        similarityMetricChar=func2str(similarityMetric);
    end
    if isequal(computing,'gpu')
        try
            canUseGPU = parallel.gpu.GPUDevice.isAvailable;
        catch ME
            canUseGPU = false;
        end
        if canUseGPU==0
            error('Please check your GPU device driver.')
        else
            g=gpuDevice();
            reset(g);
        end
        if isa(similarityMetric,'function_handle')
            similarityMetricChar=func2str(similarityMetric);
        end
        if ismember(similarityMetricChar,{'TfunctionDist','spearman'})
            warning('cannot compute distance on gpu, switching to cpu.')
            computing='cpu';
        end
    end
    [NumTFs, NumGenes] = size(RegNet);
    disp(['Learning Network with ' similarityMetricChar ' !']);
    tic;
    step = 0;
    hamming = 1;
    if isequal(precision,'single')
        TFCoop=single(TFCoop);
        RegNet=single(RegNet);
        GeneCoReg=single(GeneCoReg);
        warning off;
    end
    if isequal(computing,'gpu')
        TFCoop   = gpuArray(TFCoop);
        RegNet   = gpuArray(RegNet);
        GeneCoReg= gpuArray(GeneCoReg);
    end
    if isPUMA==1
        TFCoopInit=TFCoop;
    end
    while hamming > 0.001
        if isequal(similarityMetric,'Tfunction')
            R = Tfunction(TFCoop, RegNet);
            A = Tfunction(RegNet, GeneCoReg);
        else
            if ~isequal(similarityMetric,'minkowski')
                A = pdist2(GeneCoReg,RegNet,similarityMetric);
                R = pdist2(TFCoop, RegNet',similarityMetric);
            else
                R = pdist2(TFCoop, RegNet',similarityMetric,3);
                A = pdist2(GeneCoReg,RegNet,similarityMetric,3);
            end
            R = convertToSimilarity(R,similarityMetric);
            A = convertToSimilarity(A',similarityMetric);
        end
        
        A = respWeight*R + (1-respWeight)*A;
        if saveMemory==1
            clear R;prevDiag=diag(GeneCoReg);
            GeneCoReg=diagsquareform(GeneCoReg);
        end

        hamming = mean(abs(RegNet(:) - A(:)));
        RegNet = (1 - alpha) * RegNet + alpha * A;

        if hamming > 0.001
            if isequal(similarityMetric,'Tfunction')
                A = Tfunction(RegNet);
            else
                if ~isequal(similarityMetric,'minkowski')
                    A = pdist(RegNet,similarityMetric);
                else
                    A = pdist(RegNet,similarityMetric,3);
                end
                A = squareform(A);
                A = convertToSimilarity(A,similarityMetric);
            end
            
            A = UpdateDiagonal(A, NumTFs, alpha, step);
            TFCoop = (1 - alpha) * TFCoop + alpha * A;%clear A;
 
            if isPUMA==1
                TFCoopDiag=diag(TFCoop);
                TFCoop(sub2ind([NumTFs, NumTFs],s1,s2))=...
                    TFCoopInit(sub2ind([NumTFs, NumTFs],s1,s2));
                TFCoop(sub2ind([NumTFs, NumTFs],t1,t2))=...
                    TFCoopInit(sub2ind([NumTFs, NumTFs],t1,t2)); 
                    TFCoop(1:(size(TFCoop,1)+1):end) =TFCoopDiag; 
            end
            if isequal(similarityMetric,'Tfunction')
                A = Tfunction(RegNet');
            else
                if ~isequal(similarityMetric,'minkowski')
                    A = pdist(RegNet',similarityMetric);
                else
                    A = pdist(RegNet',similarityMetric,3);
                end
                try %sometimes MATLAB throws an out of memory error if it does have memory so we reset the gpu
                    A = squareform(A);
                catch ME
                    fprintf('caught memory error, trying to fix ... \n')
                    %save variables to disk
                    A        = gather(A);
                    GeneCoReg= gather(GeneCoReg);
                    prevDiag = gather(prevDiag);
                    TFCoop   = gather(TFCoop);
                    RegNet   = gather(RegNet);
                    hamming   = gather(hamming);
                    %reset device
                    gpuDevice(1);
                    %load variables
                    A         = gpuArray(A);
                    GeneCoReg = gpuArray(GeneCoReg);
                    prevDiag  = gpuArray(prevDiag);
                    TFCoop    = gpuArray(TFCoop);
                    RegNet    = gpuArray(RegNet);
                    hamming   = gpuArray(hamming);
                    %continue
                    A = squareform(A);
                end
                A = convertToSimilarity(A,similarityMetric);
            end
            A = UpdateDiagonal(A, NumGenes, alpha, step);
            if saveMemory==1
                % update GeneCoReg
                stdDiag = diag(A);
                A = diagsquareform(A);
                GeneCoReg = (1 - alpha) * GeneCoReg + alpha * A;
                clear A;
                GeneCoReg = squareform(GeneCoReg);
                GeneCoReg(1:(NumGenes+1):end) = alpha * stdDiag + (1 - alpha) * prevDiag;
                clear prevDiag; clear stdDiag;
            else
                GeneCoReg = (1 - alpha) * GeneCoReg + alpha * A;
            end
        end
        if verbose==1
            disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        end
        step = step + 1;
    end
    runtime = toc;
    if isequal(computing,'gpu') && gpuSave==0
        RegNet=gather(RegNet);
        gpuDevice(1);%Clear GPU device memory 
    end
    fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
    
end
