function W = otter(W,P,C)
%parameters
lambda = 0.0035;
gamma = 0.335;
Imax = 300;
eta = 0.00001;
%ADAM parameters
b1 = 0.9;
b2 = 0.999;
eps = 0.00000001;
%initial transformation
C = C/trace(C);
P = P/trace(P) + 0.0013;
W = P*W;
W = W/sqrt(trace(W*W'));

[nTF, nGenes] = size(W);
m = zeros(nTF, nGenes);
v = m;
b1t = b1;
b2t = b2;
%To save computations:
P = -P*(1-lambda) + gamma*eye(nTF);
C = -C*lambda;
%ADAM gradient descent                  
    for i = 1:Imax
       %gradient
       grad = (P*W + W*C + W*W'*W);
       m = b1*m + ((1-b1)*4)*grad;
       v = b2*v + ((1-b2)*16)*grad.^2; 
       b1t = b1t*b1;
       b2t = b2t*b2;
       alpha = sqrt(1-b2t)/(1-b1t)*eta;
       epst = eps*sqrt((1-b2t));
       %update of gene ragulatory matrix
       W = W - alpha*(m./(epst+sqrt(v)));
    end
end