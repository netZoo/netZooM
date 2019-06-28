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
% pairsNet=readtable(networkPair,'FileType','text'); % Matlab format

networkPair
pwd
fid = fopen(networkPair, 'r'); % Octave compatible format
frewind(fid);
pairsNet = textscan(fid, '%s %s %s %s', 'delimiter', '\t'); % tiny speed-up by not checking for comments
fclose(fid);

% Find number of TFs
nTFs = length(pairsNet{4}(2:end))/nGenes;

% Test if dimensions are correct
if mod(length(pairsNet{4}(2:end)),nTFs)~=0
    error('Check the number of genes!')
end

% Build network in matrix format
if prior==0
    matNet=reshape(pairsNet{4}(2:end),[nTFs,nGenes])';
elseif prior==1
    matNet=reshape(pairsNet{3}(2:end),[nTFs,nGenes])';
end

end