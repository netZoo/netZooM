function AgNet=panda_run(lib_path, exp_file, motif_file, ppi_file, panda_out, save_temp, alpha, save_pairs, modeProcess,...
               respWeight, absCoex, similarityMetric, computing, precision, verbose, saveGPUmemory)
% Description:
%               Using PANDA to infer gene regulatory network. 
%               1. Reading in input data (expression data, motif prior, TF PPI data)
%               2. Computing coexpression network
%               3. Normalizing networks
%               4. Running PANDA algorithm
%               5. Writing out PANDA network (optional)
%
% Inputs:
%               exp_file  : path to file containing gene expression as a matrix of size (g,n)
%               motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%               ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%               panda_out : path to save output PANDA network
%                           '*.txt': the final network will be saved in .txt format
%                           '*.tsv': the final network will be saved in .tsv format
%                           '*.*'  : the final network will be saved in .mat v6 format
%                           ''     : the final network will not be saved
%               save_temp : path to save updated ppi, co-expression, and gene regulation network
%                           '': the networks will not be saved
%               alpha     : learning parameter for the PANDA algorithm
%               save_pairs: (Optional) boolean parameter
%                           1:  the final network will be saved .pairs format where each line has a TF-gene edge (Cytoscape compatible)
%                           0:  the final network will not be saved in .pairs format
%               modeProcess: Refers to the procedure to filter input data.
%                           'legacy' old deprecated behavior of netZooM <0.4.1
%                                    aligns genes on gene expression and
%                                    TFs on motif.
%                           (default)'union' fills missing genes and TFs with zero rows
%                           'intersection' removes missing genes and TFs
%               respWeight:  real number between 0 and 1. Weight of the responsability matrix (default: 0.5)
%               absCoex   :  0: take the signed correlation matrix (default)
%                            1: take the absolute value of the coexpression matrix
%               similarityMetric: string containing the similarity metric that PANDA uses to 
%                           find agreement between networks. Similarity
%                           scores are kept as is, and distance scores were
%                           converted to similarities through s=1/(1+d)
%                           'Tfunction'   : (default) Modified tanimoto
%                                          similarity as described in doi:10.1371/journal.pone.0064832
%                           @TfunctionDist: same as 'Tfunction' but works
%                                          with pdist2, slower.
%                           'euclidean'  : 1/(1+ euclidean distance) 
%                           'squaredeuclidean': 1/(1+ squared euclidean distance) 
%                           'seuclidean' : 1/(1+ standardized euclidean distance) 
%                           'cityblock'  : 1/(1+ cityblock distance)
%                           'minkowski'  : 1/(1+ minkowski distance). 
%                                          with p=3. p=1,p=2, and p=inf are
%                                          covered by cityblock, euclidean,
%                                          and Chebychev respectively.
%                           'chebychev'  : 1/(1+ chebychev distance)
%                           'cosine'     : cosine of the included angle
%                           'correlation': sample correlation between points
%                           'hamming'    : 1/(1+ hamming distance)
%                           'jaccard'    : Jaccard coefficient
%                           'spearman'   : sample Spearman's rank correlation
%               computing: 'cpu'(default)
%                          'gpu' uses GPU to compute distances
%               precision: computing precision
%                          double: double precision(default)
%                          single: single precision
%               verbose  : 1 prints iterations (default)
%                          0 does not print iterations
%               saveGPUmemory: invoked only when computing='gpu'
%                              0 uses the standard implementation (default) 
%                              1 saves memory on the GPU (slower)
%                          
% Outputs:
%               AgNet     : Predicted TF-gene complete regulatory network using PANDA as a matrix of size (t,g).
%
% Authors: 
%               cychen, marieke, kglass
% 
% Notes:
%               Script adapted from Marieke's pLIONESSpcss.m, modified to run PANDA only.
% 
% Publications:
%               https://doi.org/10.1371/journal.pone.0064832 

disp(datestr(now));
% Set default parameters
if nargin <16
   saveGPUmemory=0; 
end
if nargin < 15
    verbose=1;
end
if nargin < 14
    precision='double';
end
if nargin < 13
   computing='cpu';
end
if nargin < 12
   similarityMetric='Tfunction';
end
if nargin < 11
    absCoex=0;
end
if nargin < 10
    respWeight=0.5;
end
if nargin < 9
    modeProcess='union';
end
if nargin < 8
	save_pairs=0;
end
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
[Exp,RegNet,TFCoop,TFNames,GeneNames]=processData(exp_file,motif_file,ppi_file,modeProcess);

%% ============================================================================
%% Run PANDA
%% ============================================================================
disp('Computing coexpression network:');
tic; GeneCoReg = Coexpression(Exp); toc;
if absCoex==1
    GeneCoReg=abs(GeneCoReg);
end

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
AgNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight, similarityMetric,...
    computing, precision, verbose, saveGPUmemory);

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
    if save_pairs==1
        SavePairs(TFNames, GeneNames, AgNet, RegNet, panda_out);
    end
end

disp('All done!');
disp(datestr(now));

end
