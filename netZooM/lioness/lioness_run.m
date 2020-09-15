function lioness_run(exp_file, motif_file, ppi_file, panda_file, save_dir,...
    START, END, alpha, ascii_out, lib_path, computing, saveGPUmemory, verbose, ncores)
% Description:
%             Using LIONESS to infer single-sample gene regulatory networks.
%             1. Reading in PANDA network and preprocessed middle data
%             2. Computing coexpression network
%             3. Normalizing coexpression network
%             4. Running PANDA algorithm
%             5. Writing out LIONESS networks
%             lioness_run can be called directly, otherwise for batch calls
%             1. Run PANDA first to get preprocessed middle files and aggregated PANDA network.
%             2. Set up LIONESS run-time parameters by editing `lioness_config.m`.
%             3. Run LIONESS main program via `lioness_run.sh`.
% Inputs:
%             exp_file  : path to file containing gene expression as a
%                         matrix of size (g,g) from the save_temp folder as obtained from panda_run() 
%             motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%                         from the save_temp folder as obtained from panda_run() 
%             ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%                         from the save_temp folder as obtained from panda_run() 
%             panda_file: path to the PANDA generated gene regulatory network
%             save_dir  : path to save directory. If it does not exist, it will be created.
%             START     : index of first sample to generate predicted gene regulatory network.
%             END       : index of last sample to generate predicted gene regulatory network. There will be END-START+1 single network samples generated.
%                         -1: use the index of the final gene expression sample
%             alpha     : learning parameter for the PANDA algorithm
%             ascii_out : 1 : save LIONESS networks in .txt format
%                         0 : save LIONESS networks in .mat -v6 format
%             lib_path  : path to library
%             computing: 'cpu'(default)
%                        'gpu' uses GPU to compute distances
%             saveGPUmemory: invoked only when computing='gpu'
%                            0 uses the standard implementation (default) 
%                            1 saves memory on the GPU (slower)
%             verbose  : 1 prints iterations (default)
%                        0 does not print iterations
%             ncores    : integer to specify the number of cores used to
%                         run LIONESS in parallel on ncores CPUs.
% Outputs:
%             PredNet  : Predicted single sample network as a matrix of size (t,g)
%                        This output is directly saved as a file and not
%                        as worksapce variable
% Author(s):
%             cychen, marieke, kglass
% Publications:
%             https://doi.org/10.1016/j.isci.2019.03.021

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
    fprintf('ASCII output: %d\n', ascii_out);
    addpath(lib_path);
    if nargin<14
        ncores=1;
    end
    if nargin <13
       verbose=0; 
    end
    if nargin <12
       saveGPUmemory=0; 
    end
    if nargin<11
        computing='cpu';
    end
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
    % Create the output folder if not exists
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % Sample indexes to iterate
    if END == -1
        indexes = START:NumConditions;
    else
        indexes = START:END;
    end

    if isequal(computing,'gpu')
        parpool(gpuDeviceCount); 
        parfor i = indexes
            ExpGPU = gpuArray(Exp);
            fprintf('Running LIONESS for sample %d:\n', i);
            idx = [1:(i-1), (i+1):NumConditions];  % all samples except i

            disp('Computing coexpression network:');
            tic; GeneCoReg = Coexpression(ExpGPU(idx,:)); toc;

            disp('Normalizing Networks:');
            tic; GeneCoReg = NormalizeNetwork(GeneCoReg); toc;
            GeneCoReg = gather(GeneCoReg); % it is faster to gather here although memory issues could occur with large networks
                                           % in which case, gather should after coexpression

            disp('Running PANDA algorithm:');
            LocNet = gpuPANDA(RegNet, GeneCoReg, TFCoop, alpha, 0.5,...
                'Tfunction','gpu','single',verbose, saveGPUmemory)
            PredNet = NumConditions * (AgNet - LocNet) + LocNet;

            saveGPU(PredNet,ascii_out,save_dir,i);
        end
    else
        if ncores==1
            for i = indexes
                fprintf('Running LIONESS for sample %d:\n', i);
                idx = [1:(i-1), (i+1):NumConditions];  % all samples except i

                disp('Computing coexpression network:');
                tic; GeneCoReg = Coexpression(Exp(idx,:)); toc;

                disp('Normalizing Networks:');
                tic; GeneCoReg = NormalizeNetwork(GeneCoReg); toc;

                disp('Running PANDA algorithm:');
                LocNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);
                PredNet = NumConditions * (AgNet - LocNet) + LocNet;

                disp('Saving LIONESS network:');
                if ascii_out == 1
                    f = fullfile(save_dir, sprintf('lioness.%d.txt', i));
                    tic; save(f, 'PredNet', '-ascii'); toc;
                else
                    f = fullfile(save_dir, sprintf('lioness.%d.mat', i));
                    tic; save(f, 'PredNet', '-v6'); toc;
                end
                fprintf('Network saved to %s\n', f);

                clear idx GeneCoReg LocNet PredNet f; % clean up for next run
            end
        else
            parpool(ncores);
            parfor i = indexes
                fprintf('Running LIONESS for sample %d:\n', i);
                idx = [1:(i-1), (i+1):NumConditions];  % all samples except i

                disp('Computing coexpression network:');
                tic; GeneCoReg = Coexpression(Exp(idx,:)); toc;

                disp('Normalizing Networks:');
                tic; GeneCoReg = NormalizeNetwork(GeneCoReg); toc;

                disp('Running PANDA algorithm:');
                LocNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha);
                PredNet = NumConditions * (AgNet - LocNet) + LocNet;

                saveCPU(PredNet,ascii_out,save_dir,i);
            end
        end
    end

    disp('All done!');
    disp(datestr(now));
    
end

function saveGPU(PredNet,ascii_out,save_dir,i)
% Description:
%             saveGPU circumvents the MATLAB lock on saving files inside a
%             parfor loop
% Inputs:
%             PredNet  : Predicted single sample network as a matrix of size (t,g)
%                          This output is directly saved as a file and not
%                          as worksapce variable
%             ascii_out : 1 : save LIONESS networks in .txt format
%                         0 : save LIONESS networks in .mat -v6 format
%             save_dir  : path to save directory. If it does not exist, it will be created.
%             i         : index of the network to save

    disp('Saving LIONESS network:');
    if ascii_out == 1
        f = fullfile(save_dir, sprintf('lioness.%d.txt', i));
        tic; save(f, 'PredNet', '-ascii'); toc;
    else
        f = fullfile(save_dir, sprintf('lioness.%d.mat', i));
        tic; save(f, 'PredNet', '-v6'); toc;
    end
    fprintf('Network saved to %s\n', f);
    
end

function saveCPU(PredNet,ascii_out,save_dir,i)
% Description:
%             saveCPU circumvents the MATLAB lock on saving files inside a
%             parfor loop
% Inputs:
%             PredNet  : Predicted single sample network as a matrix of size (t,g)
%                          This output is directly saved as a file and not
%                          as worksapce variable
%             ascii_out : 1 : save LIONESS networks in .txt format
%                         0 : save LIONESS networks in .mat -v6 format
%             save_dir  : path to save directory. If it does not exist, it will be created.
%             i         : index of the network to save

    disp('Saving LIONESS network:');
    if ascii_out == 1
    	f = fullfile(save_dir, sprintf('lioness.%d.txt', i));
    	tic; save(f, 'PredNet', '-ascii'); toc;
    else
    	f = fullfile(save_dir, sprintf('lioness.%d.mat', i));
    	tic; save(f, 'PredNet', '-v6'); toc;
    end
    fprintf('Network saved to %s\n', f);
    clear idx GeneCoReg LocNet PredNet f; % clean up for next run
                
end