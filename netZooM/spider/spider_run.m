function SpiderNet=spider_run(lib_path, bedtoolspath, alpha, motifhitfile,  annofile, chrinfo, ranges, regfile, outtag,motifdir, epifile,save_temp,save_pairs,spider_out,nTF )

% Description:
% 		Using SPIDER to infer epigenetically-informed gene regulatory network. 

%		  Optional steps that do not need to run everytime. The codes can be run once and parameters values can be saved
%               1. Create regulatory regions (user-defined ranges) using DefineRegulatoryRegions.m, Mostly have to do only once
%               2. Create epigenetically-filtered motif locations using CreateEpigeneticMotif.m
%               3. Create input prior network (motif prior) using BuildSPIDERprior.m
%               4. Normalize the networks
                
%		  Actual SPIDER algorithm if the SPIDER PRIOR is already constructed
  
%		            5. Running PANDA algorithm
%               6. Writing out SPIDER network (optional)


% 	Inputs for DefineRegulatoryRegions.m
%               bedtoolspath : path of the bedtools (can be installed from : "https://bedtools.readthedocs.io/en/latest/content/installation.html")
%               annofile     : file with gene annotations e.g., './ReferenceData/refseq_hg19_05292018'
%               chrinfo      : file with chromosome information e.g.,  '/InputData/ReferenceData/GenomeWideRanges.bed'
%               ranges       : user-defined input for ranges around TSS for constructing proximal or distal SPIDER networks
%   Output:
%               regfile      : path to save bedfile of regulatory regions used asinput to build SPIDER prior network. 
%	            	motifhitfile : path to file containing epigenetically informed motif information, can be created using CreateEpigeneticMotif.m



%   Input for Building SPIDER prior using BuildSPIDERprior.m 
% 
%               motifhitfile : path to file containing epigenetically informed motif information, can be created using CreateEpigeneticMotif.m
%               regfile      : path to file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m
%   Output:  
%               Adj          : path to epigenetically-filtered motif prior regulatory network for given cell line (DNase-seq data), can be created using BuildSPIDERprior.m
%               TFNames      : names of TFs in the prior network obtained from BuildSPIDERprior.m, 
%               GeneNames    : names of Genes in the prior network obtained from BuildSPIDERprior.m

%   Inputs for message-passing PANDA step in SPIDER:

%	              alpha        : parameter that determines the level of message-passing
%         		  outtag       : path to save SPIDER networks
%	            	PriorNet     : path to epigenetically-filtered motif prior regulatory network for the 
%				                      given cell line (DNase-seq data), can be created using BuildSPIDERprior.m

%   Optional PANDA Parameters for saving 
%               exp_file  : path to file containing gene expression as a matrix of size (g,g)
%               motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%               ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%               spider_out: path to save output SPIDER network
%                           '*.txt': the final network will be saved in .txt format
%                           '*.tsv': the final network will be saved in .tsv format
%                           '*.*'  : the final network will be saved in .mat v6 format
%                           ''     : the final network will not be saved
%               save_temp : path to save updated ppi, co-expression, and gene regulation network
%                           '': the networks will not be saved
%               save_pairs: (Optional) boolean parameter
%                           1:  the final network will be saved .pairs format where each line has a TF-gene edge (Cytoscape compatible)
%                           0:  the final network will not be saved in .pairs format
% 
% FINAL Output:
%               SpiderNet     : Predicted TF-gene gene complete regulatory network using SPIDER as a matrix of size (t,g).
%

% Authors: 
%               Abhijeet Sonawane, Kimberly Glass
% 
% Publications:
% 
%                
%% ============================================================================
%% Set Program Parameters and Path
%% ============================================================================

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

alpha=0.1; % level of message-passing



%%%%% PREPROCESSSING
%DefineRegulatoryRegions(annofile, ranges, regfile, chrinfo);

CreateEpigeneticMotif(epifile, motifdir, motifhitfile, bedtoolspath,nTF);

%%%% Run SPIDER %%%%

% Build SPIDER prior

[PriorNet, TFNames, GeneNames]=BuildSPIDERprior(motifhitfile, regfile, bedtoolspath);
% Run message-passing
%addpath('/udd/spaso/netZooM/netZooM/tools/')
SpiderNet=SPIDER(PriorNet, eye(length(GeneNames)), eye(length(TFNames)), alpha);

%%%% Print data to file %%%%

% reshape network information into vectors 
% TF=repmat(TFNames, 1, length(GeneNames)); TF=TF(:);
% gene=repmat(GeneNames', length(TFNames), 1); gene=gene(:);
% PriorNet=PriorNet(:);
% SpiderNet=SpiderNet(:);
% 
% % print to file
% fid=fopen([outtag, '_FinalNetwork.pairs'], 'wt');
% fprintf(fid, 'TF\tgene\tMotif\tSPIDER-prediction\n');
% for(cnt=1:length(TF))
% 	fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, PriorNet(cnt), SpiderNet(cnt));
% end
% fclose(fid);

%% ============================================================================
%% Saving SPIDER network output
%% ============================================================================
if ~isempty(spider_out)
    disp('Saving PANDA network!');
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
