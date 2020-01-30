function RegNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric,...
                    computing, precision, verbose)
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
%                           'hamming'    : 1/(1+ hamming distance) used for discrete 
%                                          variables. GPU enabled
%                           'jaccard'    : Jaccard coefficient used for discrete 
%                                          variables. GPU enabled
%                           'spearman'   : sample Spearman's rank correlation
%                           'computing'  : 'cpu'(default)
%                                          'gpu' uses GPU to compute distances
%               precision: computing precision
%                         double: double precision(default)
%                         single: single precision
%               verbose  : 1 prints iterations (Default)
%                          0 does not print iterations
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
    if nargin<8
        precision='double';
    end
    if nargin<9
        verbose=1;
    end
    if iscategorical(similarityMetric)
        similarityMetric=char(similarityMetric(1));
    end
    similarityMetricChar=similarityMetric;
    if isa(similarityMetric,'function_handle')
        similarityMetricChar=func2str(similarityMetric);
    end
    if isequal(computing,'gpu')
        RegNet = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight,...
        similarityMetric,computing,precision);
        return
    end
    if isequal(precision,'single')
        TFCoop=single(TFCoop);
        RegNet=single(RegNet);
        GeneCoReg=single(GeneCoReg);
        warning off;
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
        if verbose==1
            disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        end
        step = step + 1;
        clear A;  % release memory for next step
    end
    runtime = toc;
    fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
end