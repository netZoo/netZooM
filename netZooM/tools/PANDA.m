function RegNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric)
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
%                           'euclidean'  : 1/(1+ euclidean distance) 
%                           'squaredeuclidean': 1/(1+ squared euclidean distance) 
%                           'seuclidean' : 1/(1+ standardized euclidean distance) 
%                           'cityblock'  : 1/(1+ cityblock distance)
%                           'minkowski'  : 1/(1+ minkowski distance). 
%                                          with p=3. p=1,p=2, and p=inf are
%                                          covered by cityblock, euclidean,
%                                          and Chebychev respectively.
%                           'chebychev'  : 1/(1+ chebychev distance)
%                           'cosine'     : cosine of the included angle
%                           'correlation': sample correlation between points
%                           'hamming'    : 1/(1+ hamming distance)
%                           'jaccard'    : Jaccard coefficient
%                           'spearman'   : sample Spearman's rank correlation
%
% Outputs:
%               RegNet   : inferred gene-TF regulatory network
%
% Authors:
%               Kimberley Glass
%
% Publications:
%               https://doi.org/10.1371/journal.pone.0064832 
    if nargin<11
        respWeight=0.5;
    end
    if nargin<12
        similarityMetric='Tfunction';
    end
    [NumTFs, NumGenes] = size(RegNet);
    disp('Learning Network!');
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

function tdist = TfunctionDist(X,Y)
% Description:
%              Computes a modified version of the tanimoto distance as described
%              in Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." 
%              PloS one 8.5 (2013). https://doi.org/10.1371/journal.pone.0064832
%
% Inputs:
%               X   : 1-by-n vector
%               Y   : m2-by-n matrix
%
% Outputs:
%               tdist: m2-by-1 vector of distances
%
% Authors:
%               Kimberly Glass, cychen
%
% Notes:
%               bsxfun is more memory-efficient and faster than repmat implementation for large arrays.
%               MATLAB uses BLAS routines to do matrix multiplication. MATLAB parser recognizes X*X' as
%               a symmetric matrix multiply and will call the symmetric matrix multiply routine (only 
%               calculates about 1/2 the answer and then fills in the rest with copies, which is faster).

    tdist = X * Y';
    Bvec = sum(Y' .^ 2, 1);
    Cvec = sum(X .^ 2, 2);
    tdist = tdist ./ sqrt(bsxfun(@plus, Bvec, Cvec) - abs(tdist));

end
