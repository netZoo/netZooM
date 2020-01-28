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
    if isequal(computing,'gpu')
        try
            canUseGPU = parallel.gpu.GPUDevice.isAvailable;
        catch ME
            canUseGPU = false;
        end
        if canUseGPU==0
            error('Please check your GPU device driver.')
        end
        if ~ischar(similarityMetric)
            similarityMetricChar=func2str(similarityMetric);
        end
        if ismember(similarityMetricChar,{'TfunctionDist','spearman'})
            warning('cannot compute distance on gpu, switching to cpu.')
            computing='cpu';
        end
    end
    [NumTFs, NumGenes] = size(RegNet);
    disp('Learning Network!');
    tic;
    step = 0;
    hamming = 1;
    if isequal(computing,'gpu')
        TFCoop   = gpuArray(TFCoop);
        RegNet   = gpuArray(RegNet);
        GeneCoReg= gpuArray(GeneCoReg);
    end
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
        W = respWeight*R + (1-respWeight)*A;
        hamming = mean(abs(RegNet(:) - W(:)));
        RegNet = (1 - alpha) * RegNet + alpha * W;

        if hamming > 0.001
            if isequal(similarityMetric,'Tfunction')
                PPI = Tfunction(RegNet);
            else
                if ~isequal(similarityMetric,'minkowski')
                    PPI = pdist(RegNet,similarityMetric);
                else
                    PPI = pdist(RegNet,similarityMetric,3);
                end
                PPI = convertToSimilarity(PPI,similarityMetric);
                PPI = squareform(PPI);
            end
            PPI = UpdateDiagonal(PPI, NumTFs, alpha, step);
            TFCoop = (1 - alpha) * TFCoop + alpha * PPI;

            if isequal(similarityMetric,'Tfunction')
                CoReg2 = Tfunction(RegNet');
            else
                if ~isequal(similarityMetric,'minkowski')
                    CoReg2 = pdist(RegNet',similarityMetric);
                else
                    CoReg2 = pdist(RegNet',similarityMetric,3);
                end
                CoReg2 = convertToSimilarity(CoReg2,similarityMetric);
                CoReg2 = squareform(CoReg2);
            end
            CoReg2 = UpdateDiagonal(CoReg2, NumGenes, alpha, step);
            GeneCoReg = (1 - alpha) * GeneCoReg + alpha * CoReg2;
        end

        disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        step = step + 1;
        clear R A W PPI CoReg2;  % release memory for next step
    end
    runtime = toc;
    if isequal(computing,'gpu')
        RegNet=gather(RegNet);
    end
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
    if ~ischar(method)
        method=func2str(method);
    end
    if ismember(method,similarityList)
        mat=1-mat;
    elseif ~isequal(method,'TfunctionDist')
        mat=mat./(1+mat);
    end
end
