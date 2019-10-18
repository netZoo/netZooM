function A=DegreeAdjust(A);

k1=sum(A,1)/size(A,1);
k2=sum(A,2)/size(A,2);
A=A.*sqrt(repmat(k2,1,size(A,2)).^2+repmat(k1,size(A,1),1).^2);
