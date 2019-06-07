function normMat = NormalizeNetwork(X)
% Description:
%              Normalize an input network (matrix) X
%
% Inputs:
%               X: network adjacency matrix
%
% Outputs:
%               normMat: nromalized adjacency matrix
%
% Auhtors:
%               Kimberley Glass

    Z1 = zscore(X, 1, 1);
    Z2 = zscore(X, 1, 2);
    normMat = (Z1 + Z2) / sqrt(2);
end
