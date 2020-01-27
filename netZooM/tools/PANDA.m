function RegNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight)
% Description:
%              PANDA infers a gene regulatory network from gene expression
%              data, motif prior, and PPI between transcription factors
%
% Inputs:
%               RegNet    : motif prior of gene-TF regulatory network
%               GeneCoReg : gene-gene co-regulatory network
%               TFCoop    : PPI binding between transcription factors
%               respWeight: real number between 0 and 1. Weight of the responsability matrix (default: 0.5)
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
    [NumTFs, NumGenes] = size(RegNet);
    disp('Learning Network!');
    tic;
    step = 0;
    hamming = 1;
    while hamming > 0.001
        R = pdist2(TFCoop, RegNet',@Tfunction);
        A = pdist2(GeneCoReg,RegNet,@Tfunction);
        W = respWeight*R + (1-respWeight)*A';
        hamming = mean(abs(RegNet(:) - W(:)));
        RegNet = (1 - alpha) * RegNet + alpha * W;

        if hamming > 0.001
            PPI = pdist(RegNet,@Tfunction);
            PPI = squareform(PPI);
            PPI = UpdateDiagonal(PPI, NumTFs, alpha, step);
            TFCoop = (1 - alpha) * TFCoop + alpha * PPI;

            CoReg2 = pdist(RegNet',@Tfunction);
            CoReg2 = squareform(CoReg2);
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
