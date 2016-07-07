function normMat=NormalizeNetwork(X);

mu0=mean(X(:));
std0=std(X(:));
mu1=mean(X);
std1=std(X,1);
mu2=mean(X,2);
std2=std(X,1,2);

Z1=(X-repmat(mu1, size(X,1), 1))./repmat(std1, size(X,1), 1);
Z2=(X-repmat(mu2, 1, size(X,2)))./repmat(std2, 1, size(X,2));
normMat=Z1/sqrt(2)+Z2/sqrt(2);

%{
% checks and defaults for missing data
filter=isnan(Z2);
normMat(filter)=Z1(filter);

filter=isnan(Z1);
normMat(filter)=Z2(filter);

filter=isnan(Z1) & isnan(Z2);
normMat(filter)=(X(filter)-mu0)/std0;
%}
