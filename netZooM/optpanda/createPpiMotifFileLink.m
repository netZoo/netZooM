function [motif_file,ppi_file,pandaData]=createPpiMotifFileLink(exp_file,motifWeight,motifCutOff,...
    addCorr,motif_file,absCoex,ppi_file,thresh,oldMotif,incCoverage,...
    qpval,bridgingProteins,addChip,ctrl,ppiExp,explore)
% Description:
%               Generate a new tf-gene regulation prior for optPANDA. 
%               Please refer to the LiveScript tutorial 
%               'Accurate reconstruction of gene regulatory networks using optPANDA'
%               for a full working example. https://github.com/netZoo/netZooM/tree/master/tutorials
%
% Inputs:
%               exp_file    : path to file containing gene expression as a matrix of size (g,n)
%               motifWeight : [0,1] when addCorr==1. Weight of TF-Gene coexpression when added to the TF-gene regulation prior.
%               motifCutOff : [0,1] cutoff values in the TF-Gene coexpression, when addCorr==1
%               addCorr     : {0,1} augment the TF-gene regulation prior by adding TF-gene coexpression
%               motif_file  : TF-gene regulation prior containing p-values of FIMO runs of TF PWMs. The threshold
%                             thresh will be used to assign interactions.
%               absCoex     : {0,1} takes the absolute value of the coexpression matrix
%               ppi_file    : link to file containing TF-TF PPI
%               thresh      : [0,1] p-value threshold to assign binding in TF-gene regulation prior
%               oldMotif    : {0,1} Use a TF-gene regulation prior that
%                                   was generated previously in Lopez-Ramos
%                                   et al. https://doi.org/10.1186/s12864-017-4111-x
%               incCoverage : {0,1} increases TF coverage by adding rows of zeros in TF-gene 
%                                   regulation prior if a TF has an entry in the PPI matrix
%               qpval       : {0,1} use p-value or corrected p-value for TF-gene binding motif
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
%               explore     : 0 saves the optimized motif and PPI files
%                             1 does not save result file, meant for
%                               exploring the solution space
%
% Outputs:
%               motif_file  : tf-gene regulation prior for optPANDA in
%                             .pairs format.
%               ppi_file    : tf-tf ppi prior for optPANDA in
%                             .pairs format.
%
% Authors: 
%               Marouen Ben Guebila 01/2020

    % read expression and generate coexpression network
    disp('Creating an optimal motif and PPI priors!');
    expTbl=readtable(exp_file,'FileType','text');
    Exp = expTbl{:,2:end};
    GeneNames = expTbl{:,1};
    [NumGenes, NumConditions] = size(Exp);
    fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
    Exp = Exp';  % transpose expression matrix from gene-by-sample to sample-by-gene
        
    % read ppi file to get TF names
    ppitbl      = readtable(ppi_file,'FileType','text');
    TFNames= unique(ppitbl.Var1);
    NumTFs = length(TFNames);
    TFCoop = zeros(NumTFs);
    TF1    = ppitbl.Var1;
    TF2    = ppitbl.Var2;
    weight = ppitbl.Var3;
    [~,i]  = ismember(TF1, TFNames);
    [~,j]  = ismember(TF2, TFNames);
    TFCoop(sub2ind([NumTFs, NumTFs], i, j)) = weight;
    TFCoop(sub2ind([NumTFs, NumTFs], j, i)) = weight; 

    %read precomputed results
    if isstruct(motif_file)
        if isequal(motif_file.selectedMotif,'Hugo_motifCellLine.txt')
            RegNettmp    = motif_file.a.RegNet;
            tmpTFNames   = motif_file.a.TFNames;
            tmpGeneNames = motif_file.a.GeneNames;
        elseif isequal(motif_file.selectedMotif,'regMatPval1e3.txt')
            RegNettmp    = motif_file.b.RegNet;
            tmpTFNames   = motif_file.b.TFNames;
            tmpGeneNames = motif_file.b.GeneNames;
        elseif isequal(motif_file.selectedMotif,'regMatQval005.txt')
            RegNettmp    = motif_file.c.RegNet;
            tmpTFNames   = motif_file.c.TFNames;
            tmpGeneNames = motif_file.c.GeneNames;
        elseif isequal(motif_file.selectedMotif,'motif_complete_reduced.txt')
            RegNettmp    = motif_file.d.RegNet;
            tmpTFNames   = motif_file.d.TFNames;
            tmpGeneNames = motif_file.d.GeneNames;
        end
        motif_file=motif_file.selectedMotif;
        [ii,i]   = ismember(tmpTFNames, TFNames);
        [jj,j]   = ismember(tmpGeneNames, GeneNames);
        i        = i(i~=0);ii=find(ii);
        j        = j(j~=0);jj=find(jj);
        RegNet   = zeros(NumTFs, NumGenes);
        RegNet(i,j) = RegNettmp(ii,jj);
    else
        [TF, gene, weight] = textread(motif_file, '%s%s%f');
        [~,i]   = ismember(TF, TFNames);
        [~,j]   = ismember(gene, GeneNames);
        RegNet  = zeros(NumTFs, NumGenes);
        indComm = i&j;
        i       = i(indComm);
        j       = j(indComm);
        weight  = weight(indComm);
        RegNet(sub2ind([NumTFs, NumGenes], i, j)) = weight;
        fprintf('%d TFs and %d edges!\n', NumTFs, length(weight));
    end
    
    % Apply threshold for p values only, q-values are thresholded at 0.05
    if oldMotif==0 && incCoverage==1 && qpval==0        
        RegNet(RegNet <= thresh & RegNet > 0)  = 1;
        RegNet(RegNet > thresh  & RegNet < 1 ) = 0;
    end
    
    % add chip-seq data
    if addChip==1
        [gt]=readGt('CHIP_SEQ_REMAP_ALL.txt');
        [~,igt,itf]=intersect(gt{:,1},TFNames);
        RegNet(itf,:)=0;
        for i=1:length(igt)
            [~,~,ic]  = intersect(gt{igt(i),2:end}, GeneNames);
            RegNet(itf(i),ic)= 1;
        end
    elseif addChip==2
        [gt]=readGt('CHIP_SEQ_REMAP_ALL.txt');
        [~,igt,itf]=intersect(gt{:,1},TFNames);
        RegNet(itf,:)=-1;
        for i=1:length(igt)
            [~,~,ic]  = intersect(gt{igt(i),2:end}, GeneNames);
            RegNet(itf(i),ic)= 2;
        end
    end
    
    % compute gene-gene coexpression
    GeneCoReg = Coexpression(Exp);
    
    % compute new ppi by scaling ppi with tf-tf coexpression
    if ppiExp==1
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = zeros(NumTFs,NumTFs);
        multMotif(ipr,ipr) = abs(GeneCoReg(ifi,ifi));
        TFCoop             = TFCoop.*multMotif;
    elseif ppiExp==2
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = ones(NumTFs,NumTFs);
        multMotif(ipr,ipr) = abs(GeneCoReg(ifi,ifi));
        TFCoop             = TFCoop.*multMotif;  
    elseif ppiExp==3
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = zeros(NumTFs,NumTFs);
        multMotif(ipr,ipr) = abs(GeneCoReg(ifi,ifi));
        multMotif(multMotif == 0) = mean2(multMotif);
        TFCoop             = TFCoop.*multMotif; 
    elseif ppiExp==4
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = zeros(NumTFs,NumTFs);
        multMotif(ipr,ipr) = GeneCoReg(ifi,ifi);
        TFCoop             = TFCoop.*multMotif;
    elseif ppiExp==5
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = ones(NumTFs,NumTFs);
        multMotif(ipr,ipr) = GeneCoReg(ifi,ifi);
        TFCoop             = TFCoop.*multMotif;  
    elseif ppiExp==6
        [~,ifi,ipr]   = intersect(GeneNames, TFNames);
        multMotif          = zeros(NumTFs,NumTFs);
        multMotif(ipr,ipr) = GeneCoReg(ifi,ifi);
        multMotif(multMotif == 0) = mean2(multMotif);
        TFCoop             = TFCoop.*multMotif; 
    end
    
    % add tf-gene coexpression to motif prior
    if addCorr==1
        % Find TF index in gene names
        addMotif=zeros(NumTFs,NumGenes);
        % Remove diagonal to avoid bias for self loops
        GeneCoReg=GeneCoReg-diag(diag(GeneCoReg));
        for i=1:NumTFs
            [setInt,ifi,~]=intersect(GeneNames, TFNames(i));
            if ~isempty(setInt)
                addMotif(i,:)=GeneCoReg(ifi,:);
            end
        end
        addMotif = abs(addMotif); % we take the absolute value of the correlation matrix
        addMotif(addMotif == 0) = mean2(addMotif);
        % apply parameters
        addMotif(addMotif < motifCutOff) = 0;
        RegNet = (1-motifWeight)*RegNet+motifWeight*abs(addMotif);
    elseif addCorr==2
        % Find TF index in gene names
        addMotif=zeros(NumTFs,NumGenes);
        for i=1:NumTFs
            [setInt,ifi,~]=intersect(GeneNames, TFNames(i));
            if ~isempty(setInt)
                addMotif(i,:)=GeneCoReg(ifi,:);
            end
        end
        addMotif = abs(addMotif); % we take the absolute value of the correlation matrix
        addMotif(addMotif == 0) = mean2(addMotif);
         % apply parameters
        addMotif(addMotif < motifCutOff) = 0;
        RegNet = (1-motifWeight)*RegNet+motifWeight*abs(addMotif);
    elseif addCorr==3
        % Find TF index in gene names
        addMotif=zeros(NumTFs,NumGenes);
        % Remove diagonal to avoid bias for self loops
        GeneCoReg=GeneCoReg-diag(diag(GeneCoReg));
        for i=1:NumTFs
            [setInt,ifi,~]=intersect(GeneNames, TFNames(i));
            if ~isempty(setInt)
                addMotif(i,:)=GeneCoReg(ifi,:);
            end
        end
        addMotif = abs(addMotif); % we take the absolute value of the correlation matrix
         % apply parameters
        addMotif(addMotif < motifCutOff) = 0;
        RegNet = (1-motifWeight)*RegNet+motifWeight*abs(addMotif);
    elseif addCorr==4
        % Find TF index in gene names
        addMotif=zeros(NumTFs,NumGenes);
        for i=1:NumTFs
            [setInt,ifi,~]=intersect(GeneNames, TFNames(i));
            if ~isempty(setInt)
                addMotif(i,:)=GeneCoReg(ifi,:);
            end
        end
        addMotif = abs(addMotif); % we take the absolute value of the correlation matrix
        % apply parameters
        addMotif(addMotif < motifCutOff) = 0;
        RegNet = (1-motifWeight)*RegNet+motifWeight*abs(addMotif);
    end
    
    if explore==1
        %keep results in memory
        pandaData.RegNet=RegNet;pandaData.GeneCoReg=GeneCoReg;
        pandaData.TFCoop=TFCoop;
        ppi_file='';
        motif_file='';
    elseif explore==0%save result files and give link
        display('Saving optimized PPI and Motif priors')
        % Create motif file name
        [~,name,~]=fileparts(exp_file);
        motif_file = [name '_' motif_file(1:end-4) '_MW' num2str(motifWeight) '_MC' num2str(motifCutOff)...
            '_AC' num2str(addCorr) '_ABS' num2str(absCoex) '_THR' num2str(thresh)...
            '_OM' num2str(oldMotif) '_IC' num2str(incCoverage) '_QP' num2str(qpval)...
            '_BR' num2str(bridgingProteins) '_CHIP' num2str(addChip) ...
            '_PE' num2str(ctrl) '.txt'];
        [~,name,~]=fileparts(ppi_file);
        ppi_file = [name '_' ppi_file(1:end-4) '_MW' num2str(motifWeight) '_MC' num2str(motifCutOff)...
            '_AC' num2str(addCorr) '_ABS' num2str(absCoex) '_THR' num2str(thresh)...
            '_OM' num2str(oldMotif) '_IC' num2str(incCoverage) '_QP' num2str(qpval)...
            '_BR' num2str(bridgingProteins) '_CHIP' num2str(addChip) ...
            '_PE' num2str(ctrl) '.txt'];

        % Save motif file
        TF    = repmat(TFNames, 1, length(GeneNames));
        gene  = repmat(GeneNames', length(TFNames), 1);
        TF    = TF(:);
        gene  = gene(:);
        RegNet= RegNet(:);%linearizes by column
        % Saving file
        fid   = fopen(motif_file, 'wt');
        for cnt=1:length(TF)
            fprintf(fid, '%s\t%s\t%f\n', TF{cnt}, gene{cnt}, RegNet(cnt) );
        end
        fclose(fid);

        % Save PPI file
        TF    = repmat(TFNames, 1, length(TFNames));
        TF2   = repmat(TFNames', length(TFNames), 1);
        TF    = TF(:);
        TF2   = TF2(:);
        TFCoop= TFCoop(:);%linearizes by column
        % Saving file
        fid   = fopen(ppi_file, 'wt');
        for cnt=1:length(TF)
            fprintf(fid, '%s\t%s\t%f\n', TF{cnt}, TF2{cnt}, TFCoop(cnt) );
        end
        fclose(fid);
    end
end
