function tdist = TfunctionDist(X,Y)
% Description:
%             Computes a modified version of the tanimoto distance as described
%             in Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." 
%             PloS one 8.5 (2013). https://doi.org/10.1371/journal.pone.0064832
% Inputs:
%             X   : 1-by-n vector
%             Y   : m2-by-n matrix
% Outputs:
%             tdist: m2-by-1 vector of distances
% Author(s):
%             Kimberly Glass, cychen
% Notes:
%             bsxfun is more memory-efficient and faster than repmat implementation for large arrays.
%             MATLAB uses BLAS routines to do matrix multiplication. MATLAB parser recognizes X*X' as
%             a symmetric matrix multiply and will call the symmetric matrix multiply routine (only 
%             calculates about 1/2 the answer and then fills in the rest with copies, which is faster).

    tdist = X * Y';
    Bvec = sum(Y' .^ 2, 1);
    Cvec = sum(X .^ 2, 2);
    tdist = tdist ./ sqrt(bsxfun(@plus, Bvec, Cvec) - abs(tdist));

end