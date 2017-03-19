%% Using LIONESS to infer single-sample gene regulatory networks.
%%
%% 1. Reading in PANDA network and preprocessed middle data
%% 2. Computing coexpression network
%% 3. Normalizing coexpression network
%% 4. Running PANDA algorithm
%% 5. Writing out LIONESS networks
%%
%% Authors: cychen, marieke, kglass

disp(datestr(now));

%% ============================================================================
%% Set Program Parameters and Path
%% ============================================================================
% Run configuration to set parameter first (e.g., run('lioness_config.m');)
fprintf('Input expression file: %s\n', exp_file);
fprintf('Input motif file: %s\n', motif_file);
fprintf('Input PPI file: %s\n', ppi_file);
fprintf('Input PANDA network: %s\n', panda_file);
fprintf('Output LIONESS folder: %s\n', save_dir);
fprintf('Sample index: %d - %d\n', START, END);
fprintf('Alpha: %.2f\n', alpha);
addpath(lib_path);

%% ============================================================================
%% Read in Data
%% ============================================================================
disp('Reading in expression data!');
tic
    X = load(exp_file);
    Exp = X.Exp;
    [NumConditions, NumGenes] = size(Exp);  % transposed expression
    fprintf('%d genes and %d conditions!\n', NumGenes, NumConditions);
toc

disp('Reading in motif data!');
tic
    X = load(motif_file);
    RegNet = X.RegNet;
toc

disp('Reading in ppi data!');
tic
    X = load(ppi_file);
    TFCoop = X.TFCoop;
toc
 
disp('Reading in PANDA network!');
tic
    X = load(panda_file);
    AgNet = X.AgNet;
toc

%% ============================================================================
%% Run LIONESS
%% ============================================================================
% Create output folder if not exists
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for i = START:END
    fprintf('Running LIONESS for sample %d:\n', i);
    idx = [1:(i-1), (i+1):NumConditions];  % all samples except i

    disp('Computing coexpresison network:');
    tic; GeneCoReg = Coexpression(Exp(idx,:)); toc;

    disp('Normalizing Networks:');
    tic; GeneCoReg = NormalizeNetwork(GeneCoReg); toc;

    disp('Running PANDA algorithm:');
    LocNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);
    PredNet = NumConditions * (AgNet - LocNet) + LocNet;

    disp('Saving LIONESS network:');
    f = fullfile(save_dir, sprintf('lioness.%d.mat', i));
    tic; save(f, 'PredNet', '-v6'); toc;
    fprintf('Network saved to %s\n', f);

    % Ascii output for testing
    %f = fullfile(save_dir, sprintf('lioness.%d.txt', i));
    %x = PredNet(:);
    %tic; save(f, 'x', '-ascii'); toc;
    %fprintf('Network saved to %s\n', f);

    % For low memory machine
    %clear idx GeneCoReg LocNet PredNet f; % clean up for next run
end

disp('All done!');
disp(datestr(now));
