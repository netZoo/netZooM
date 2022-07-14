function [GeneMotif,GeneNamesExp,TfMotif,TFNamesInit,NumConditions,ExpInit,TF,gene,weightMotif,weightPPI,TF1,TF2,SampleNames]=readData(exp_file,motif_file,ppi_file)
% Description:
%             readData reads the input files for PANDA.
% Inputs:
%             exp_file  : file for gene expression
%             motif_file: file for motif data
%             ppi_file  : file for TF PPI data
% Ouputs:
%             GeneMotif    : Gene names from gene motif data
%             GeneNamesExp : Gene names from gene expression data
%             TfMotif      : TF names from gene motif data
%             TFNamesInit  : TF names from gene PPI data
%             NumConditions: Number of gene expression data
%             ExpInit      : Gene expression matrix
%             TF           : list of TF edges in motif (source)
%             gene         : list of gene edges in motif (target)
%             weightMotif  : edge weight in the motif network
%             weightPPI    : edge weight in the PPI network
%             TF1          : list of TF edges in PPI (source)
%             TF2          : list of TF edges in PPI (target)
%             SampleNames  : gene expression samples IDs
% Author:     
%             Marouen Ben Guebila 12/2019

    % Read expression
    disp('Reading in expression data!');
    tic
        exp_file_tbl = readtable(exp_file,'FileType','text','PreserveVariableNames',1);
        SampleNames = exp_file_tbl.Properties.VariableNames;
        ExpInit      = exp_file_tbl{:,2:end};
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