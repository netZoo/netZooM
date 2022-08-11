function [TF,gene,AgNetCol,RegNetCol]=SavePairs(TFNames, GeneNames, AgNet, RegNet, outtag, header)
% Description:
%             A function to save a complete graph in matrix format to a pairs format where each line represents an edge in the network.
%             The output file will have as much lines as edges in the network and will have the predicted edge weights as well as the 
%             binary edge weights from the prior motif data.
% Inputs:
%             TFNames  : names of t TFs
%             GeneNames: names of g genes
%             AgNet    : predicted gene regulation network using PANDA of size (t,g)
%             RegNet   : prior gene regulation network obtained using TF motif scan of size (t,g)
%                        []: saves only AgNet edges
%             outtag   : name of saved file 
%                        '': does not save file on disk
%             header   : (0/1) write header for the outtag
%
% Outputs:
%             TF       : Column list of TFs (source)
%             gene     : Column list of genes (target)
%             AgNetCol : Edge weight in predicted gene regulatory network
%             RegNetCol: Edge weight in prior gene regulatory network
%
% Author(s):
%             Kimberly Glass, Marouen Ben Guebila

    if nargin<6
        header=0;
    end
    % Check input dimensions (Results in misannotated files otherwise)
    if size(GeneNames,1)==1
        GeneNames=GeneNames';
    end
    % Reshape network information into vectors and print to file
    TF    = repmat(TFNames, 1, length(GeneNames));
    gene  = repmat(GeneNames', length(TFNames), 1);
    TF    = TF(:);
    gene  = gene(:);    
    AgNetCol = AgNet(:);
    if isempty(RegNet)==0
        RegNetCol= RegNet(:);
    else
        RegNetCol=[];
    end

    % Save file
    if strcmp(outtag,'')==0
        fid   = fopen([outtag, '_FinalNetwork.pairs'], 'wt');
        if isempty(RegNet)==0
            if header==1
                fprintf(fid, 'TF\tgene\tMotif\tPANDA-prediction\n');
            end
            for cnt=1:length(TF)
                fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, RegNetCol(cnt), AgNetCol(cnt));
            end
        else
            if header==1
                fprintf(fid, 'TF\tgene\tPANDA-prediction\n');
            end
            for cnt=1:length(TF)
                fprintf(fid, '%s\t%s\t%f\n', TF{cnt}, gene{cnt}, AgNetCol(cnt));
            end
        end
        fclose(fid);
    end

end
