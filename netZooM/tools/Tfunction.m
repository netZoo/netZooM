function Amat = Tfunction(X,Y)
% Description:
%              Updates Amat matrix using matrices X and Y         
%
% Inputs:
%               X   : adjacency matrix
%               Y   : adjacency matrix
%
% Outputs:
%               Amat: updated matrix from X and Y
%
% Authors:
%               Kimberly Glass, cychen
%
% Notes:
%               bsxfun is more memory-efficient and faster than repmat implementation for large arrays.
%               MATLAB uses BLAS routines to do matrix multiplication. MATLAB parser recognizes X*X' as
%               a symmetric matrix multiply and will call the symmetric matrix multiply routine (only 
%               calculates about 1/2 the answer and then fills in the rest with copies, which is faster).
    switch nargin
        case 1
            Cvec = bsxfun(@plus, sum(X .^ 2, 2)', sum(X .^ 2, 2));
            Amat = X * X';
            Amat = Amat ./ sqrt(Cvec - abs(Amat));
        case 2
            Amat = X * Y;
            Bvec = sum(Y .^ 2, 1);
            Cvec = sum(X .^ 2, 2);
            Amat = Amat ./ sqrt(bsxfun(@plus, Bvec, Cvec) - abs(Amat));
    end
end
