function matNet=Pairs2Mat(networkPair,nGenes,prior)
% Paris2Mat transforms a TF-Gene network in .pairs network to a complete
% matrix. It can also save the prior in matrix format.
%
% Inputs:
%         networkPair: path to network in .pairs format with nGenes*nTFs
%                      edges
%         nGenes     : number of genes in network
%         prior      : 0: matNet is a matrix of the final network
%                      1: matNet is a matrix of the prior network
%                   
% Outputs:
%         matNet     : network in nGenes by nTfs matrix format
%
% Author: Marouen Ben Guebila 6/19

% Read network in .pairs format
pairsNet=readtable(networkPair,'FileType','text');

% Find number of TFs
nTFs = length(pairsNet.Var4)/nGenes;

% Test if dimensions are correct
if mod(pairsNet.Var4,nTFs)~=0
    error('Check the number of genes!')
end

% Build network in matrix format
if prior==0
    matNet=reshape(pairsNet.Var4,[nTFs,nGenes])';
elseif prior==1
    matNet=reshape(pairsNet.Var3,[nTFs,nGenes])';
end

end