% Normalize an input network (matrix) X
function normMat = NormalizeNetwork(X)
    Z1 = zscore(X, 1, 1);
    Z2 = zscore(X, 1, 2);
    normMat = (Z1 + Z2) / sqrt(2);
end
