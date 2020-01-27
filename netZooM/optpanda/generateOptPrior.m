function [motif_file,GeneNames,allTFName,ppi_file]=generateOptPrior(exp_file,incCoverage,bridgingProteins,oldMotif,...
    qpval,precomputed,motif_fil,motifWeight,motifCutOff,addCorr,absCoex,thresh,...
    addChip,ctrl,ppiExp)
% Description:
%               Generate a new tf-gene regulation prior for optPANDA. 
%
% Inputs:
%               exp_file    : path to file containing gene expression as a matrix of size (g,n)
%               incCoverage : {0,1} increases TF coverage by adding rows of zeros in TF-gene 
%                                   regulation prior if a TF has an entry in the PPI matrix
%               bridgingProteins : computation of higher order PPI interaction network through linking 
%                                  TFs if they sahre a non-TF neighbor.
%                                  PPIs were taken from STRINGdb v 11.0
%                                  0 : PPI with all STRINGdb links (binding, catalysis, etc ...)
%                                  1 : 0 + addition of edge if two TFs
%                                      share a neighbor. This will account for
%                                      first order neighbors or bridging
%                                      proteins in TF complexes.
%                                  2 : 1 + addition of edge it two TFs
%                                      share neighors who share neighbors i.e.,
%                                      second order bridging proteins
%                                  3 : PPI with only edges that account for
%                                      binding interactions
%                                  4 : 3 + first-order neighbors
%                                  5 : 4 + second-order neighbors
%                                  6 : 5 + third-order neighbors
%                                  7 : 6 + fourth-order neighbors
%               oldMotif    : {0,1} Use a TF-gene regulation prior that
%                                   was generated previously in Lopez-Ramos
%                                   et al. https://doi.org/10.1186/s12864-017-4111-x
%               qpval       : {0,1} use p-value or corrected p-value for TF-gene binding motif
%               precomputed : {0,1} use presaved results
%               motif_fil   : {0,1} use p-value or corrected p-value for TF-gene binding motif
%               motifWeight : [0,1] when addCorr==1. Weight of TF-Gene coexpression when added to the TF-gene regulation prior.
%               motifCutOff : [0,1] cutoff values in the TF-Gene coexpression, when addCorr==1
%               addCorr     : {0,1} augment the TF-gene regulation prior by adding TF-gene coexpression
%               absCoex     : {0,1} takes the absolute value of the coexpression matrix
%               thresh      : [0,1] p-value threshold to assign binding in TF-gene regulation prior
%               addChip     : {0,1,2} adds chip-seq in the TF-gene regulation prior
%                                  1 : adds chipseq data from all of remap (union of all chipseq experiments
%                                      i.e., nonspecific to conetxt) 
%                                  2 : same as 1 but uses 2 and -1 for
%                                      presence/absence to differentiate from
%                                      PWMs that use 1 and 0.
%               ctrl        : {0,1,2,3} Dummy variable used in the optimisation process 
%                                       to compute the performance of a null variable 
%               ppiExp      : {0,1,2,3,4,5,6} scaling of PPI data by TF-TF
%                              Coexpression (tfco)
%                             1 : PPI*.|tfco| + fill missing with 0
%                             2 : PPI*.|tfco| + fill missing with 1
%                             3 : PPI*.|tfco| + fill missing with mean
%                             4 : PPI*.tfco + fill missing with 0
%                             5 : PPI*.tfco + fill missing with 1
%                             6 : PPI*.tfco + fill missing with mean
%                          
% Outputs:
%               motif_file : tf-gene regulation prior for optPANDA in
%                            .pairs format.
%               GeneNames  : list of genes in motif_file
%               allTFName  : list of TF names in motif_file
%
% Authors: 
%               Marouen Ben Guebila 01/2020

    % fetch gene names
    a   = readtable(exp_file,'FileType','text');
    Exp = a{:,2:end};
    GeneNames = a{:,1};
     
    if incCoverage==1
        switch bridgingProteins
            %/!\ Important ppi with all links has a diagonal of 1_nTfs.
            % ppi with binding only has a diagonal of 0_nTFs.
            case 0
                ppi_file = 'ppi_complete.txt';
            case 1
                ppi_file = '1ppi_complete.txt';
            case 2
                ppi_file = '2ppi_complete.txt';
            case 3
                ppi_file = 'ppi_complete_bind.txt';
            case 4
                ppi_file = '1ppi_complete_bind.txt';
            case 5
                ppi_file = '2ppi_complete_bind.txt';
            case 6
                ppi_file = '3ppi_complete_bind.txt';
            case 7
                ppi_file = '4ppi_complete_bind.txt';
        end
        [TF1, TF2, weight] = textread(ppi_file, '%s%s%f');
        % read TF names
        allTFName  = unique(TF1);
        disp('Reading in motif data with extra coverage!');
        if oldMotif==1
            motif_file = 'motif_complete_reduced.txt';
        elseif oldMotif==0
            if qpval==1 % use FDR-corrected p-value
                motif_file = 'regMatQval005.txt';
            elseif qpval==0 % use p-value
                motif_file = 'regMatPval1e3.txt';
            end
        end
    else
        motif_file = 'Hugo_motifCellLine.txt';
        ppi_file   = 'ppi2015_freezeCellLine.txt';
        [TF1, TF2, weight] = textread(ppi_file, '%s%s%f');
        allTFName  = unique(TF1);
    end
    if precomputed == 0 
       if isstruct(motif_fil)
           motif_file=motif_fil; 
       end
    end
    if precomputed == 0
        % create motif file
        [motif_file,ppi_file]=createPpiMotifFileLink(exp_file,motifWeight,motifCutOff,...
            addCorr,motif_file,absCoex,ppi_file,thresh,oldMotif,...
            incCoverage,qpval,bridgingProteins,addChip,ctrl,ppiExp);
    elseif precomputed==1
        % provide link to motif file based on parameters
        [filepath,name,ext]=fileparts(exp_file);
        motif_file = [name '_' motif_file(1:end-4) '_MW' num2str(motifWeight) '_MC' num2str(motifCutOff)...
            '_AC' num2str(addCorr) '_ABS' num2str(absCoex) '_THR' num2str(thresh)...
            '_OM' num2str(oldMotif) '_IC' num2str(incCoverage) '_QP' num2str(qpval)...
            '_BR' num2str(bridgingProteins) '_CHIP' num2str(addChip) ...
            '_PE' num2str(ctrl) '.txt'];
        [filepath,name,ext]=fileparts(ppi_file);
        ppi_file = [name '_' ppi_file(1:end-4) '_MW' num2str(motifWeight) '_MC' num2str(motifCutOff)...
        '_AC' num2str(addCorr) '_ABS' num2str(absCoex) '_THR' num2str(thresh)...
        '_OM' num2str(oldMotif) '_IC' num2str(incCoverage) '_QP' num2str(qpval)...
        '_BR' num2str(bridgingProteins) '_CHIP' num2str(addChip) ...
        '_PE' num2str(ctrl) '.txt'];
    end
    
end