function RegNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric,...
    computing)
% Description:
%              PANDA infers a gene regulatory network from gene expression
%              data, motif prior, and PPI between transcription factors
%
% Inputs:
%               RegNet    : motif prior of gene-TF regulatory network
%               GeneCoReg : gene-gene co-regulatory network
%               TFCoop    : PPI binding between transcription factors
%               respWeight: real number between 0 and 1. Weight of the responsability matrix (default: 0.5)
%               similarityMetric: string containing the similarity metric that PANDA uses to 
%                           find agreement between networks. Similarity
%                           scores are kept as is, and distance scores were
%                           converted to similarities through s=1/(1+d)
%                           'Tfunction'   : (Default) Modified tanimoto
%                                          similarity as described in doi:10.1371/journal.pone.0064832
%                           @TfunctionDist: same as 'Tfunction' but works
%                                          with pdist2, slower.
%                           'euclidean'  : 1/(1+ euclidean distance). GPU enabled
%                           'squaredeuclidean': 1/(1+ squared euclidean
%                                               distance). GPU enabled
%                           'seuclidean' : 1/(1+ standardized euclidean
%                                          distance) . GPU enabled
%                           'cityblock'  : 1/(1+ cityblock distance). GPU enabled
%                           'minkowski'  : 1/(1+ minkowski distance). 
%                                          with p=3. p=1,p=2, and p=inf are
%                                          covered by cityblock, euclidean,
%                                          and Chebychev respectively. GPU enabled
%                           'chebychev'  : 1/(1+ chebychev distance). GPU enabled
%                           'cosine'     : cosine of the included angle.GPU enabled
%                           'correlation': sample correlation between
%                                          points. GPU enabled.
%                           'hamming'    : 1/(1+ hamming distance). GPU enabled
%                           'jaccard'    : Jaccard coefficient. GPU enabled
%                           'spearman'   : sample Spearman's rank correlation
%                           'computing'  : 'cpu'(default)
%                                          'gpu' uses GPU to compute distances
%
% Outputs:
%               RegNet   : inferred gene-TF regulatory network
%
% Authors:
%               Kimberley Glass
%
% Publications:
%               https://doi.org/10.1371/journal.pone.0064832 
    if nargin<5
        respWeight=0.5;
    end
    if nargin<6
        similarityMetric='Tfunction';
    end
    if nargin<7
        computing='cpu';
    end
    if iscategorical(similarityMetric)
        similarityMetric=char(similarityMetric(1));
    end
    similarityMetricChar=similarityMetric;
    if isequal(computing,'gpu')
        RegNet = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight,...
        similarityMetric,computing);
        return
    end
    [NumTFs, NumGenes] = size(RegNet);
    disp(['Learning Network with ' similarityMetricChar ' !']);
    tic;
    step = 0;
    hamming = 1;
    while hamming > 0.001
        if isequal(similarityMetric,'Tfunction')
            R = Tfunction(TFCoop, RegNet);
            A = Tfunction(RegNet, GeneCoReg);
        else
            if ~isequal(similarityMetric,'minkowski')
                R = pdist2(TFCoop, RegNet',similarityMetric);
                A = pdist2(GeneCoReg,RegNet,similarityMetric);
            else
                R = pdist2(TFCoop, RegNet',similarityMetric,3);
                A = pdist2(GeneCoReg,RegNet,similarityMetric,3);
            end
            R = convertToSimilarity(R,similarityMetric);
            A = convertToSimilarity(A',similarityMetric);
        end
        A = respWeight*R + (1-respWeight)*A;
        clear R;

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
                A = convertToSimilarity(A,similarityMetric);
                A = squareform(A);
            end
            A = UpdateDiagonal(A, NumTFs, alpha, step);
            TFCoop = (1 - alpha) * TFCoop + alpha * A;

            if isequal(similarityMetric,'Tfunction')
                A = Tfunction(RegNet');
            else
                if ~isequal(similarityMetric,'minkowski')
                    A = pdist(RegNet',similarityMetric);
                else
                    A = pdist(RegNet',similarityMetric,3);
                end
                A = convertToSimilarity(A,similarityMetric);
                A = squareform(A);
            end
            A = UpdateDiagonal(A, NumGenes, alpha, step);
            GeneCoReg = (1 - alpha) * GeneCoReg + alpha * A;
        end

        disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        step = step + 1;
        clear A;  % release memory for next step
    end
    runtime = toc;
    fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
end

function mat=convertToSimilarity(mat,method)
% Description:
%              convertToSimilarity converts a distance to similarity
%
% Inputs:
%               mat   : n-by-m distance matrix
%               method: the metric being converted
%                       'cosine','correlation','jaccard','spearman': will
%                       be 1-distance.
%                       'euclidean','seuclidean','squaredeuclidean','cityblock'
%                       'minkowski','chebychev','hamming': will be 1./(mat+1)
%
% Outputs:
%               mat : n-by-m similarity matrix

    similarityList={'cosine','correlation','jaccard','spearman'};
    if isa(method,'function_handle')
        method=func2str(method);
    end
    if ismember(method,similarityList)
        mat=1-mat;
    elseif ~isequal(method,'TfunctionDist')
        mat=mat./(1+mat);
    end
end

function RegNet = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric,...
    computing)
% Description:
%              GPU-accelerated PANDA. 
%
% Inputs:
%               RegNet    : motif prior of gene-TF regulatory network
%               GeneCoReg : gene-gene co-regulatory network
%               TFCoop    : PPI binding between transcription factors
%               respWeight: real number between 0 and 1. Weight of the responsability matrix (default: 0.5)
%               similarityMetric: string containing the similarity metric that PANDA uses to 
%                           find agreement between networks. Similarity
%                           scores are kept as is, and distance scores were
%                           converted to similarities through s=1/(1+d)
%                           'Tfunction'   : (Default) Modified tanimoto
%                                          similarity as described in doi:10.1371/journal.pone.0064832
%                           @TfunctionDist: same as 'Tfunction' but works
%                                          with pdist2, slower.
%                           'euclidean'  : 1/(1+ euclidean distance). GPU enabled
%                           'squaredeuclidean': 1/(1+ squared euclidean
%                                               distance). GPU enabled
%                           'seuclidean' : 1/(1+ standardized euclidean
%                                          distance) . GPU enabled
%                           'cityblock'  : 1/(1+ cityblock distance). GPU enabled
%                           'minkowski'  : 1/(1+ minkowski distance). 
%                                          with p=3. p=1,p=2, and p=inf are
%                                          covered by cityblock, euclidean,
%                                          and Chebychev respectively. GPU enabled
%                           'chebychev'  : 1/(1+ chebychev distance). GPU enabled
%                           'cosine'     : cosine of the included angle.GPU enabled
%                           'correlation': sample correlation between
%                                          points. GPU enabled.
%                           'hamming'    : 1/(1+ hamming distance). GPU enabled
%                           'jaccard'    : Jaccard coefficient. GPU enabled
%                           'spearman'   : sample Spearman's rank correlation
%                           'computing'  : 'cpu'(default)
%                                          'gpu' uses GPU to compute distances
%
% Outputs:
%               RegNet   : inferred gene-TF regulatory network
%
% Authors:
%               Kimberley Glass
%
% Publications:
%               https://doi.org/10.1371/journal.pone.0064832 
    if nargin<5
        respWeight=0.5;
    end
    if nargin<6
        similarityMetric='Tfunction';
    end
    if nargin<7
        computing='cpu';
    end
    if iscategorical(similarityMetric)
        similarityMetric=char(similarityMetric(1));
    end
    similarityMetricChar=similarityMetric;
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
    if isequal(computing,'gpu')
        TFCoop   = gpuArray(TFCoop);
        RegNet   = gpuArray(RegNet);
        GeneCoReg= gpuArray(GeneCoReg);
    end
    while hamming > 0.001
        %GeneCoReg=squareformdiag(GeneCoReg) and set diag 
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
        clear R;stdDiag = (1-alpha)* diag(GeneCoReg);
        GeneCoReg=diagsquareform(GeneCoReg);clear GeneCoReg;

        hamming = mean(abs(RegNet(:) - A(:)));
        RegNet = (1 - alpha) * RegNet + alpha * A;

        if hamming > 0.001
            if isequal(similarityMetric,'Tfunction')
                A = Tfunction(RegNet);
                A = diagsquareform(A);
            else
                if ~isequal(similarityMetric,'minkowski')
                    A = pdist(RegNet,similarityMetric);
                else
                    A = pdist(RegNet,similarityMetric,3);
                end
                A = convertToSimilarity(A,similarityMetric);
            end
            A = UpdateDiagonal(A, NumTFs, alpha, step);
            TFCoop = (1 - alpha) * TFCoop + alpha * A;

            if isequal(similarityMetric,'Tfunction')
                A = Tfunction(RegNet');
            else
                if ~isequal(similarityMetric,'minkowski')
                    A = pdist(RegNet',similarityMetric);
                else
                    A = pdist(RegNet',similarityMetric,3);
                end
                A = convertToSimilarity(A,similarityMetric);
                %A = squareform(A);
            end
            if 1
                GeneCoReg = diagsquareform(GeneCoReg);
                GeneCoReg = (1 - alpha) * GeneCoReg + alpha * A;
                clear A;
                GeneCoReg = squareform(GeneCoReg);
                GeneCoReg = UpdateDiagonal(GeneCoReg, NumGenes, alpha, step);
                GeneCoReg(1:NumGenes+1:end) = alpha * diag(GeneCoReg);
                GeneCoReg(1:NumGenes+1:end) = stdDiag + diag(GeneCoReg);
            end
        end

        disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        step = step + 1;
    end
    runtime = toc;
    if isequal(computing,'gpu')
        RegNet=gather(RegNet);
    end
    fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
end
