function matNet=Pairs2Mat(networkPair,nGenes,prior)
% Description:
%             Pairs2Mat transforms a TF-Gene network in .pairs network to a complete
%             matrix. It can also save the prior in matrix format
% Inputs:
%             networkPair: path to network in .pairs format with nGenes*nTFs
%                          edges
%             nGenes     : number of genes in network
%             prior      : 0: matNet is a matrix of the final network
%                          1: matNet is a matrix of the prior network                  
% Outputs:
%             matNet     : network in nGenes by nTfs matrix format
% Author(s):
%             Marouen Ben Guebila 6/19

    % Read network in .pairs format
    % pairsNet=readtable(networkPair,'FileType','text'); % Matlab format
    fid = fopen(networkPair, 'r'); % Octave compatible format
    frewind(fid);
    pairsNet = textscan(fid, '%s %s %s %s', 'delimiter', '\t');
    fclose(fid);

    if isempty(str2num(pairsNet{4}{1})) % check if variable name exist 
        start=2;
    else
        start=1;
    end
    % Find number of TFs
    nTFs = length(pairsNet{4}(start:end))/nGenes;

    % Test if dimensions are correct
    if mod(length(pairsNet{4}(start:end)),nTFs)~=0
        error('Check the number of genes!')
    end

    % Build network in matrix format
    if prior==0
        matNet=reshape(pairsNet{4}(start:end),[nTFs,nGenes])';
    elseif prior==1
        matNet=reshape(pairsNet{3}(start:end),[nTFs,nGenes])';
    end

    % Convert to double
    matNet=str2double(matNet);
    
end
