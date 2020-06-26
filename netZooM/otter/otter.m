function W = otter(W,P,C,lambda,gamma,Imax,eta)
% Description:
%               OTTER infers gene regulatory networks using TF DNA binding
%               motif (W), TF PPI (P), and gene coexpression (C) through 
%               minimzing the following objective:
%                                  min f(W) 
%               with f(W) = (1-lambda)*||WW' - P||^2 + lambda*||W'W - C||^2 + (gamma/2)*||W||^2
%
% Inputs:
%               W     : TF-gene regulatory network based on TF motifs as a
%                       matrix of size (t,g), g=number of genes, t=number of TFs
%               P     : TF-TF protein interaction network as a matrix of size (t,t)
%               C     : gene coexpression as a matrix of size (g,g) 
%               lambda: tuning parameter in [0,1] (higher gives more weight to C)
%               gamma : regularization parameter
%               Imax  : number of iterations of the algorithm
%               eta   : learning rate
%               bexp  : exponent influencing learning rate (higher means smaller)
%
%
% Outputs:
%               W  : Predicted TF-gene complete regulatory network as an adjacency matrix of size (t,g).
%global parameters
if nargin<4
    lambda = 0.0035;
end
if nargin<5
    gamma = 0.335;
end
if nargin<6
    Imax = 32;
end
if nargin<7
    eta = 0.00001;
end
%ADAM parameters
b1 = 0.9;
b2 = 0.999;
eps = 0.00000001;
%initial transformation
C = C/trace(C);
P = P+2.2;
P = P/trace(P);
W = P*W;
W = W/trace(W*W');

[t, g] = size(W);
m = zeros(t, g);
v = m;
b1t = b1;
b2t = b2;
%To save computations:
P = -P*(1-lambda) + gamma*eye(t);
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
    %update of gene regulatory matrix
    W = W - alpha*(m./(epst+sqrt(v)));
end
end