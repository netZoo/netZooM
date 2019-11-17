function A=DegreeAdjust(A);
% Description:
%      Function to adjust the degree so that the hub nodes are not penalized in z-score transformation
%
% Input 
%	   A : Input Network
%
% Output 
%	   A: Output Network	
 
    k1=sum(A,1)/size(A,1);
    k2=sum(A,2)/size(A,2);
    A=A.*sqrt(repmat(k2,1,size(A,2)).^2+repmat(k1,size(A,1),1).^2);
end