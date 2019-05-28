function []=panda_run(exp_file,motif_file,ppi_file,panda_out,save_temp,alpha)
% Description:
%               Using PANDA to infer gene regulatory network. 
%               1. Reading in input data (expression data, motif prior, TF PPI data)
%               2. Computing coexpression network
%               3. Normalizing networks
%               4. Running PANDA algorithm
%               5. Writing out PANDA network (optional)
%
% Inputs:
%               exp_file  : path to file containing gene expression
%               motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs
%               ppi_file  : path to file containing TF-TF interaction graph
%               panda_out : path to save output PANDA network
%                           '*.txt': the final network will be saved in .txt format
%                           '*.tsv': the final network will be saved in .tsv format
%                           '*.*'  : the final network will be saved in .mat v6 format
%                           ''     : the final network will not be saved
%               save_temp : path to save updated ppi, co-expression, and gene regulation network
%                           '': the networks will not be saved
%               alpha     : learning parameter for the PANDA algorithm
%
% 
% Authors: 
%               cychen, marieke, kglass
% 
% Notes:
%               Script adapted from Marieke's pLIONESSpcss.m, modified to run PANDA only.
% 
% Publications:
%               https://doi.org/10.1371/journal.pone.0064832
%

disp(datestr(now));

%% ============================================================================
%% Set Program Parameters and Path
%% ============================================================================
% Run configuration to set parameter first (e.g., run('panda_config.m');)
fprintf('Input expression file: %s\n', exp_file);
fprintf('Input motif file: %s\n', motif_file);
fprintf('Input PPI file: %s\n', ppi_file);
fprintf('Output PANDA network: %s\n', panda_out);
fprintf('Output temp data: %s\n', save_temp);
fprintf('Alpha: %.2f\n', alpha);
addpath(lib_path);

%% ============================================================================
%% Read in Data
%% ============================================================================
disp('Reading in expression data!');
tic
    fid = fopen(exp_file, 'r');
    headings = fgetl(fid);
    n = length(regexp(headings, '\t'));
    frewind(fid);
    %Exp = textscan(fid, ['%s', repmat('%f', 1, n)], 'delimiter', '\t', 'CommentStyle', '#');
    Exp = textscan(fid, ['%s', repmat('%f', 1, n)], 'delimiter', '\t'); % tiny speed-up by not checking for comments
    fclose(fid);
    GeneNames = Exp{1};
    Exp = cat(2, Exp{2:end});
    [NumGenes, NumConditions] = size(Exp);
    fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
    Exp = Exp';  % transpose expression matrix from gene-by-sample to sample-by-gene
toc

disp('Reading in motif data!');
tic
    [TF, gene, weight] = textread(motif_file, '%s%s%f');
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
        [TF1, TF2, weight] = textread(ppi_file, '%s%s%f');
        [~,i] = ismember(TF1, TFNames);
        [~,j] = ismember(TF2, TFNames);
        TFCoop(sub2ind([NumTFs, NumTFs], i, j)) = weight;
        TFCoop(sub2ind([NumTFs, NumTFs], j, i)) = weight;
        fprintf('%d PPIs!\n', length(weight));
    end
toc
 
% Clean up variables to release memory
clear headings n GeneNames TF gene TF1 TF2 weight;

%% ============================================================================
%% Run PANDA
%% ============================================================================
disp('Computing coexpresison network:');
tic; GeneCoReg = Coexpression(Exp); toc;

disp('Normalizing Networks:');
tic
    RegNet = NormalizeNetwork(RegNet);
    GeneCoReg = NormalizeNetwork(GeneCoReg);
    TFCoop = NormalizeNetwork(TFCoop);
toc

if ~isempty(save_temp)
    disp('Saving the transposed expression matrix and normalized networks:');
    if ~exist(save_temp, 'dir')
        mkdir(save_temp);
    end
    tic
        save(fullfile(save_temp, 'expression.transposed.mat'), 'Exp', '-v7.3');  % 2G+
        save(fullfile(save_temp, 'motif.normalized.mat'), 'RegNet', '-v6');  % fast
        save(fullfile(save_temp, 'ppi.normalized.mat'), 'TFCoop', '-v6');  % fast
    toc
end

clear Exp;  % Clean up Exp to release memory (for low-memory machine)

disp('Running PANDA algorithm:');
AgNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);

%% ============================================================================
%% Saving PANDA network output
%% ============================================================================
if ~isempty(panda_out)
    disp('Saving PANDA network!');
    tic
        [pathstr, name, ext] = fileparts(panda_out);
        switch ext
            case '.txt'
                save(panda_out, 'AgNet', '-ascii');
            case '.tsv'
                save(panda_out, 'AgNet', '-ascii', '-tabs');
            otherwise
                save(panda_out, 'AgNet', '-v6');
        end
    toc
end

disp('All done!');
disp(datestr(now));

end
