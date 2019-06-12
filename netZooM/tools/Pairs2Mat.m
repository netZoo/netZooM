function matNet=Pairs2Mat(networkPair,nGenes)
% Paris2Mat transforms a TF-Gene network in .pairs network to a complete
% matrix
%
% Inputs:
%         networkPair: path to network in .pairs format with nGenes*nTFs
%                      edges
%         nGenes     : number of genes in network
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
if ~isinteger(nTFs)
    error('Check the number of genes!')

% Build network in matrix format
matNet=reshape(pairsNet.Var4,[nTFs,nGenes])';

end