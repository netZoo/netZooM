% Implementation notes:
% bsxfun is more memory-efficient and faster than repmat implementation for large arrays.
% MATLAB uses BLAS routines to do matrix multiplication. MATLAB parser recognizes X*X' as
% a symmetric matrix multiply and will call the symmetric matrix multiply routine (only 
% calculates about 1/2 the answer and then fills in the rest with copies, which is faster).
function Amat = Tfunction(X,Y)
    switch nargin
        case 1
            Amat = X * X';
            Cvec = sum(X .^ 2, 2);
            Amat = Amat ./ sqrt(bsxfun(@plus, Cvec', Cvec) - abs(Amat));
        case 2
            Amat = X * Y;
            Bvec = sum(Y .^ 2, 1);
            Cvec = sum(X .^ 2, 2);
            Amat = Amat ./ sqrt(bsxfun(@plus, Bvec, Cvec) - abs(Amat));
end
