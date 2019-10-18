% Description:HAVE TO EDIT THE FOLLOWIN PART
%               Using SPIDER to infer gene regulatory network. 
%               1. Reading in input data (expression data, motif prior, TF PPI data)
%               2. Computing coexpression network
%               3. Normalizing networks
%               4. Running PANDA algorithm
%               5. Writing out PANDA network (optional)
%
% Inputs:
%               exp_file  : path to file containing gene expression as a matrix of size (g,g)
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
% 
% Outputs:
%               AgNet     : Predicted TF-gene gene complete regulatory network using PANDA as a matrix of size (t,g).
%
% Authors: 
%               Abhijeet Sonawane, Kimberly Glass
% 
% 
% Publications:
% SOME DAY
%               zhttps://doi.org/10.1371/journal.pone.0064832 


%%%% PARAMETER REGION %%%%

% SPIDER parameters -- note this area can be extended to include gene expression and PPI information
alpha=0.1; % level of message-passing
outtag='OutputNetworks/GM12878_ProximalSPIDER'; % name of file to print the final network info
motifhitfile='InputData/EpiMotifFiles/GM12878_filtered_motiflocations.bed'; % file storing epigenetically informed motif information, can be created with CreateEpigeneticMotif.m
regfile='InputData/RegulatoryRegions/ProximalRegulatoryRegions.bed'; % file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m

% location of codes needed to run SPIDER
bedtoolspath=''; % where bedtools is installed, set equal to '' if bedtools is already on the system path
funcpath='./SPIDER_v2/'; % path to where all the SPIDER functions are stored
addpath(funcpath);


%%%% Data Preprocessing Codes %%%%

% If needed, create the regfile with DefineRegulatoryRegions.m
annofile='InputData/ReferenceData/refseq_hg19_05292018'; % file with gene annotations
chrinfo='InputData/ReferenceData/GenomeWideRanges.bed'; % file with chromosome information
ranges={[-1000,+500]};
DefineRegulatoryRegions(annofile, ranges, regfile, chrinfo);

% If needed, create the motifhit file using CreateEpigeneticMotif.m
motifdir='InputData/MotifBedFiles/'; % where the original motif scan files are stored (one bed file per motif)
epifile='InputData/DNaseBedFiles/GM12878_DnasePeaks.bed'; % file with open chromatin regions
CreateEpigeneticMotif(epifile, motifdir, motifhitfile, bedtoolspath);

%%%% Run SPIDER %%%%

% Build SPIDER prior
[PriorNet, TFNames, GeneNames]=BuildSPIDERprior(motifhitfile, regfile, bedtoolspath);
% Run message-passing
SpiderNet=SPIDER(PriorNet, eye(length(GeneNames)), eye(length(TFNames)), alpha);

%%%% Print data to file %%%%

% reshape network information into vectors 
TF=repmat(TFNames, 1, length(GeneNames)); TF=TF(:);
gene=repmat(GeneNames', length(TFNames), 1); gene=gene(:);
PriorNet=PriorNet(:);
SpiderNet=SpiderNet(:);

% print to file
fid=fopen([outtag, '_FinalNetwork.pairs'], 'wt');
fprintf(fid, 'TF\tgene\tMotif\tSPIDER-prediction\n');
for(cnt=1:length(TF))
	fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, PriorNet(cnt), SpiderNet(cnt));
end
fclose(fid);
