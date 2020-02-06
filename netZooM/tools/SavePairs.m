function SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag)
% Description:
%             A function to save a complete graph in matrix format to a pairs format where each line represents an edge in the network.
%             The output file will have as much lines as edges in the network and will have the predicted edge weights as well as the 
%             binary edge weights from the prior motif data.
%
% Inputs:
%             TFNames  : names of t TFs
%             GeneNames: names of g genes
%             AgNet    : predicted gene regulation network using PANDA of size (t,g)
%             RegNet   : prior gene regulation network obtained using TF motif scan of size (t,g)
%             outtag   : name of saved file 
%
% Authors:
%            Kimberley Glass, Marouen Ben Guebila

% Reshape network information into vectors and print to file
    TF    = repmat(TFNames, 1, length(GeneNames));
    gene  = repmat(GeneNames', length(TFNames), 1);
    TF    = TF(:);
    gene  = gene(:);
    RegNet= RegNet(:);
    AgNet = AgNet(:);

    fid   = fopen([outtag, '_FinalNetwork.pairs'], 'wt');
    fprintf(fid, 'TF\tgene\tMotif\tPANDA-prediction\n');
    for(cnt=1:length(TF))
        fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, RegNet(cnt), AgNet(cnt));
    end
    fclose(fid);

end
