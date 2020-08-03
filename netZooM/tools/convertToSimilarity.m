function mat=convertToSimilarity(mat,method)
% Description:
%            convertToSimilarity converts a distance to similarity
% Inputs:
%             mat   : n-by-m distance matrix
%             method: the metric being converted
%                     'cosine','correlation','jaccard','spearman': will
%                     be 1-distance.
%                     'euclidean','seuclidean','squaredeuclidean','cityblock'
%                     'minkowski','chebychev','hamming': will be 1./(mat+1)
%                     Tfunction and TfunctionDist remain the same
% Outputs:
%             mat : n-by-m similarity matrix

    similarityList={'cosine','correlation','jaccard','spearman'};
    if isa(method,'function_handle')
        method=func2str(method);
    end
    if ismember(method,similarityList)
        mat=1-mat;
    elseif ~ismember(method,{'TfunctionDist','Tfunction'})
        mat=mat./(1+mat);
    end
    
end