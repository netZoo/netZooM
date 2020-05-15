function onCoex=onlineCoexpression(si,n,mi,std,cov)
% Description:
%              onlineCoexpression computes the correlation matrix of n
%              samples deprived of sample si, using the correlation matrix
%              of n samples. ~4-8x faster when number of genes ~ number of
%              observations and 12x-35x when number of samples is much
%              larger than the variables.
%              Particularly intersting in iterative computation of several
%              coexpression matrix with large samlples or when computing large 
%              matrices that require consequent GPU/CPU memory.
%
% Inputs:
%               si : k samples to remove as a k by genes vector
%               n  : number of all the samples
%               mi : mean of all the samples
%               std: std of all the samples
%               cov: covariance matrix of all the samples 
% 
% Outputs:
%               onCoex : co-expression matrix deprived of samples si
%
% Author(s):
%               Marouen Ben Guebila 5/2020
%

    [k,m]=size(si);
    % First we compute the new mean online
    newm   = (1/(n-1)) * ( (mi*n)- si);
    
    % Then we compute the new std online using the orthogonality trick
    newstd = sqrt( (std.^2 - (1/n) * (si - newm).^2 ) * (n-1)/(n-2) );
    
    % Then we compute the new covariance online
    onCov= (1/(n-2)) * ( (cov*(n-1)) - ( (n/(n-1)) * ((si-mi)'*(si-mi)) ) );
    
    % Finally, we derive the new coexpression online
    onCoex= onCov ./ (newstd'*newstd);
    
    % We set the diagonal explicitly to avoid numerical stability
    onCoex(1:m+1:end)=1; 

end