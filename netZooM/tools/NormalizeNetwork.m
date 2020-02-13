function normMat = NormalizeNetwork(X)
% Description:
%              Normalize an input network (matrix) X
%
% Inputs:
%               X: network adjacency matrix
%
% Outputs:
%               normMat: normalized adjacency matrix
%
% Auhtors:
%               Kimberley Glass

    %compute means and stds
    mu2 = mean(X,2);
    std2= std(X,1,2);
    mu1=mean(X);
    std1=std(X,1);
    mu0=mean(X(:));
    std0=std(X(:));
    normMat = (zscore(X, 1, 1) + zscore(X, 1, 2)) / sqrt(2);
    %MATLAB assigns zeros for the zscores of null std
    f2=all(zscore(X, 1, 2) == 0,2);
    f1=all(zscore(X, 1, 1) == 0,1);
    
    % checks and defaults for missing data
    Z0=(X-mu0)/std0;
    Z2=(X(:,f1)-repmat(mu2, 1, size(X(:,f1),2)))./repmat(std2, 1, size(X(:,f1),2));
    normMat(:,f1)=(Z2+Z0(:,f1))/sqrt(2);clear Z2;
    Z1=(X(f2,:)-repmat(mu1, size(X(f2,:),1), 1))./repmat(std1, size(X(f2,:),1), 1);
    normMat(f2,:)=(Z1+Z0(f2,:))/sqrt(2);clear Z1;
    linInd=(find(f1)-1)*size(normMat,1)+find(f2);
    normMat(linInd)=2*Z0(linInd)/sqrt(2);
end
