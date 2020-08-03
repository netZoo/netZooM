function diagMat = UpdateDiagonal(diagMat, num, alpha, step)
% Description:
%             Updates the diagonal of diagMat
% Inputs:
%             diagMat: diagonal matrix
%             num    : integer
%             alpha  : learning rate of the PANDA algorithm
%             step   : step number in the PANDA algorithm
% Outputs:
%             diagMat: updated diagonal matrix
% Author(s):
%             Kimberley Glass

    diagMat(1:(num+1):end) = nan;
    diagstd = nanstd(diagMat, 1);
    diagMat(1:(num+1):end) = diagstd * num * exp(2 * alpha * step);
    
end
