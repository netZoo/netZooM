% Compute gene-gene coexpression network for a sample-by-gene matrix X
% Note that each gene is a column in X
function GeneCoReg = Coexpression(X)
    GeneCoReg = corr(X, 'type', 'pearson', 'rows', 'pairwise');
    % Detecting nan in the coexpression network
    % e.g., genes with no expression variation across samples
    if any(any(isnan(GeneCoReg), 2))
        NumGenes = size(GeneCoReg, 1);
        GeneCoReg(1:NumGenes+1:NumGenes^2) = 1;  % set the diagonal to 1
        GeneCoReg(isnan(GeneCoReg)) = 0; % set nan to 0
    end
end
