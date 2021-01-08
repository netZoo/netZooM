function SpiderNet=spider_run(lib_path, bedtoolspath, alpha, motifhitfile,...  
                              annofile, chrinfo, ranges, regfile, outtag, motifdir,... 
                              epifile, save_temp, save_pairs, spider_out, nTF, computing)
% Description:
% 		Using SPIDER to infer epigenetically-informed gene regulatory network. 
%		Optional steps that do not need to run everytime. The codes can be run once and parameters values can be saved
%             1. Create regulatory regions (user-defined ranges) using DefineRegulatoryRegions.m, Mostly have to do only once
%             2. Create epigenetically-filtered motif locations using CreateEpigeneticMotif.m
%             3. Create input prior network (motif prior) using BuildSPIDERprior.m
%             4. Normalize the networks              
%		      Actual SPIDER algorithm if the SPIDER PRIOR is already constructed
%		      5. Running PANDA algorithm
%             6. Writing out SPIDER network (optional)
% Inputs: 
%             Parameters of message-passing PANDA step in SPIDER:
%             alpha        : parameter that determines the level of message-passing
%         	  outtag       : path to save SPIDER networks
%	          PriorNet     : path to epigenetically-filtered motif prior regulatory network for the 
%				                      given cell line (DNase-seq data), can be created using BuildSPIDERprior.m
%             Optional PANDA Parameters for saving: 
%             exp_file  : path to file containing gene expression as a matrix of size (g,g)
%             motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%             ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%             spider_out: path to save output SPIDER network
%                         '*.txt': the final network will be saved in .txt format
%                         '*.tsv': the final network will be saved in .tsv format
%                         '*.*'  : the final network will be saved in .mat v6 format
%                         ''     : the final network will not be saved
%             save_temp : path to save updated ppi, co-expression, and gene regulation network
%                         '': the networks will not be saved
%             save_pairs: (Optional) boolean parameter
%                         1:  the final network will be saved .pairs format where each line has a TF-gene edge (Cytoscape compatible)
%                         0:  the final network will not be saved in .pairs format
%             computing    : 'cpu'(default)
%                            'gpu' uses GPU to compute distances
% Outputs:
%             SpiderNet     : Predicted TF-gene gene complete regulatory network using SPIDER as a matrix of size (t,g).
% Author(s):
%             Abhijeet Sonawane, Kimberly Glass

    %% ============================================================================
    %% Set Program Parameters and Path
    %% ============================================================================
    if nargin < 16
       computing='cpu';
    end
    % Run configuration to set parameter first (e.g., run('panda_config.m');)

    fprintf('Input epigenetically informed motif file: %s\n', motifhitfile);
    fprintf('Input containing regulatory regions for genes: %s\n', regfile);
    fprintf('Input gene annotations file: %s\n', annofile);
    fprintf('Input chromosome information file: %s\n', chrinfo);
    fprintf('Input Ranges file: %s\n', ranges);
    fprintf('Input Original Motif Scans file: %s\n', motifdir);
    fprintf('Input PPI file: %s\n', epifile);

    % Run configuration to set parameter first (e.g., run('panda_config.m');)

    fprintf('Output SPIDER network: %s\n', spider_out);
    fprintf('Output temp data: %s\n', save_temp);
    fprintf('Alpha: %.2f\n', alpha);
    addpath(lib_path);

    %%%% PARAMETER REGION %%%%

    % SPIDER parameters -- note this area can be extended to include gene expression and PPI information

    %alpha=0.1; % level of message-passing

    %%%%% PREPROCESSSING
    if(~isempty(ranges))
        DefineRegulatoryRegions(annofile, ranges, regfile, chrinfo);
    end

    CreateEpigeneticMotif(epifile, motifdir, motifhitfile, bedtoolspath,nTF);

    %%%% Run SPIDER %%%%

    % Build SPIDER prior

    [PriorNet, TFNames, GeneNames]=BuildSPIDERprior(motifhitfile, regfile, bedtoolspath);

    %temporary reduction in number of genes for TRAVIS
    numGenes = 100 
    GeneNames = GeneNames(1:numGenes);
    PriorNet = PriorNet(:,1:numGenes);
    % Remove this section for including all genes.


    % Run message-passing
    SpiderNet=SPIDER(PriorNet, eye(length(GeneNames)), eye(length(TFNames)), alpha, computing);

    %% ============================================================================
    %% Saving SPIDER network output
    %% ============================================================================
    if ~isempty(spider_out)
        disp('Saving SPIDER network!');
        tic
            [pathstr, name, ext] = fileparts(spider_out);
            switch ext
                case '.txt'
                    save(spider_out, 'SpiderNet', '-ascii');
                case '.tsv'
                    save(spider_out, 'SpiderNet', '-ascii', '-tabs');
                otherwise
                    save(spider_out, 'SpiderNet', '-v6');
            end
        toc
        if save_pairs==1
            SavePairs(TFNames, GeneNames, SpiderNet, RegNet, spider_out);
        end
    end

    disp('All done!');
    disp(datestr(now));
    
end
