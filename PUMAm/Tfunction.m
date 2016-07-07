function Amat=Tfunction(X,Y);

Amat=(X*Y);
Bmat=repmat(sum(Y.^2,1), size(X,1), 1);
Cmat=repmat(sum(X.^2,2), 1, size(Y,2));

Amat=Amat./sqrt(Bmat+Cmat-abs(Amat));
