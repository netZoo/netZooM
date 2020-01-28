function [Exp,RegNet,TFCoop,TFNames,GeneNames]=processData(exp_file,motif_file,ppi_file,modeProcess)
% Description:
%               processData process the input data before running PANDA in
%               three different modes.
%
% Inputs:
%               exp_file  : path to file containing gene expression as a matrix of size (g,g)
%               motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%               ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%               modeProcess: 'legacy' old deprecated behavior of netZooM <= 0.5
%                            (Default)'union' fills missing genes and TFs with zero rows
%                            'intersection' removes missing genes and TFs
% Outputs:   
%               Exp       : aligned expression matrix
%               RegNet    : aligned motif prior matrix
%               TFCoop    : aligned PPI matrix
%               TFNames   : transcirption factor names
%               GeneNames : gene names
%
% Author:       Marouen Ben Guebila 12/2019

    if isequal(modeProcess,'legacy')
        disp('Reading in expression data!');
        tic
            exp_file_tbl=readtable(exp_file,'FileType','text');
            Exp = exp_file_tbl{:,2:end};
            GeneNames = exp_file_tbl{:,1};
            [NumGenes, NumConditions] = size(Exp);
            fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
            Exp = Exp';  % transpose expression matrix from gene-by-sample to sample-by-gene
        toc

        disp('Reading in motif data!');
        tic
            fid = fopen(motif_file, 'r');
            C = textscan(fid, '%s%s%f');
            fclose(fid);
            TF=C{1,1};gene=C{1,2};weight=C{1,3};
            TFNames = unique(TF);
            NumTFs = length(TFNames);
            [~,i] = ismember(TF, TFNames);
            [~,j] = ismember(gene, GeneNames);
            RegNet = zeros(NumTFs, NumGenes);
            RegNet(sub2ind([NumTFs, NumGenes], i, j)) = weight;
            fprintf('%d TFs and %d edges!\n', NumTFs, length(weight));
        toc

        disp('Reading in ppi data!');
        tic
            TFCoop = eye(NumTFs);
            if(~isempty(ppi_file))
                fid = fopen(ppi_file, 'r');
                C = textscan(fid, '%s%s%f');
                fclose(fid);
                TF1=C{1,1};TF2=C{1,2};weight=C{1,3};
                [~,i] = ismember(TF1, TFNames);
                [~,j] = ismember(TF2, TFNames);
                TFCoop(sub2ind([NumTFs, NumTFs], i, j)) = weight;
                TFCoop(sub2ind([NumTFs, NumTFs], j, i)) = weight;
                fprintf('%d PPIs!\n', length(weight));
            end
        toc
    elseif isequal(modeProcess,'union')
        [GeneMotif,GeneNamesExp,TfMotif,TFNamesInit,NumConditions,...
            ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2]=...
            readData(exp_file,motif_file,ppi_file);
        GeneNames=unique(union(GeneMotif,GeneNamesExp));
        TFNames  =unique(union(TfMotif,TFNamesInit));
        [Exp,RegNet,TFCoop]=populateData(GeneNames,TFNames,NumConditions,...
            GeneNamesExp,ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2);
    elseif isequal(modeProcess,'intersection')
        [GeneMotif,GeneNamesExp,TfMotif,TFNamesInit,NumConditions,...
            ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2]=...
            readData(exp_file,motif_file,ppi_file);
        GeneNames=intersect(GeneMotif,GeneNamesExp);
        TFNames  =intersect(TfMotif,TFNamesInit);
        [Exp,RegNet,TFCoop]=populateData(GeneNames,TFNames,NumConditions,...
            GeneNamesExp,ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2);
    end
end

function [Exp,RegNet,TFCoop]=populateData(GeneNames,TFNames,NumConditions,...
    GeneNamesExp,ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2)
    NumTFs=length(TFNames);NumGenes=length(GeneNames);
    %Initialize result
    RegNet   = zeros(NumTFs,NumGenes);
    Exp      = zeros(NumGenes,NumConditions);
    TFCoop   = zeros(NumTFs,NumTFs);
    %Populate result
    %Gene expression
    [~,ig]= ismember(GeneNamesExp,GeneNames);
    ig
    size(ExpInit)
    Exp(find(ig),:)= ExpInit(find(ig),:);
    Exp      = Exp'; 
    %Motif
    [~,i] = ismember(TF, TFNames);
    [~,j] = ismember(gene, GeneNames);
    indCommMotif = i&j;
    i      = i(indCommMotif);
    j      = j(indCommMotif);
    weightMotif = weightMotif(indCommMotif);
    RegNet(sub2ind([NumTFs, NumGenes], i, j)) = weightMotif;
    fprintf('%d TFs and %d edges!\n', NumTFs, length(weightMotif));
    %PPI
    [~,i] = ismember(TF1, TFNames);
    [~,j] = ismember(TF2, TFNames);
    indCommPPI = i&j;
    i      = i(indCommPPI);
    j      = j(indCommPPI);
    weightPPI = weightPPI(indCommPPI);
    TFCoop(sub2ind([NumTFs, NumTFs], i, j)) = weightPPI;
    TFCoop(sub2ind([NumTFs, NumTFs], j, i)) = weightPPI;    
end

function [GeneMotif,GeneNamesExp,TfMotif,TFNamesInit,NumConditions,...
            ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2]=...
            readData(exp_file,motif_file,ppi_file)
    % Read expression
    % Read expression
    disp('Reading in expression data!');
    tic
        exp_file_tbl=readtable(exp_file,'FileType','text');
        ExpInit = exp_file_tbl{:,2:end};
        GeneNamesExp = exp_file_tbl{:,1};
        [NumGenes, NumConditions] = size(ExpInit);
        fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
    
    toc
    if length(unique(GeneNamesExp)) ~= length(GeneNamesExp)
        error('There are duplicate genes in the expression matrix.')
    end
    % Read motif
    disp('Reading in motif data!');
    motif_file_id = fopen(motif_file);
    C = textscan(motif_file_id, '%s%s%f');
    fclose(motif_file_id);
    TF=C{1,1};gene=C{1,2};weightMotif=C{1,3};
    TfMotif  = unique(TF);
    GeneMotif= unique(gene);
    % Read PPI
    disp('Reading in ppi data!');
    if(~isempty(ppi_file))
        ppi_file_id = fopen(ppi_file);
        C = textscan(ppi_file_id, '%s%s%f');
        fclose(motif_file_id);
        TF1=C{1,1};TF2=C{1,2};weightPPI=C{1,3};
    end
    TFNamesInit=unique(TF1);
    if ~isequal(TFNamesInit,unique(TF2))
        error('PPI data has missing information.')
    end   
end