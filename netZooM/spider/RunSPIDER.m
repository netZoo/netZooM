% Description:
%               Using SPIDER to infer epigenetically-informed gene regulatory network. 
%               1. Create regulatory regions (user-defined ranges) using DefineRegulatoryRegions.m
%		2. Create epigenetically-filtered motif locations using CreateEpigeneticMotif.m
%		3. Create input prior network (motif prior) using BuildSPIDERprior.m
%               4. Normalizing networks
%               4. Running PANDA algorithm
%               5. Writing out SPIDER network (optional)
%
% Inputs:
%               alpha        : parameter that determines the level of message-passing
%		motifhitfile : path to file containing epigenetically informed motif information, can be created using CreateEpigeneticMotif.m
%		regfile      : path to file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m
%		PriorNet     : path to epigenetically-filtered motif prior regulatory network for given cell line (DNase-seq data), can be created using BuildSPIDERprior.m
%               TFNames      : names of TFs in the prior network obtained from BuildSPIDERprior.m, 
%		GeneNames    : names of Genes in the prior network obtained from BuildSPIDERprior.m
%		outtag       : path to save SPIDER networks	

% 
% Outputs:
%               SpiderNet     : Predicted TF-gene gene complete regulatory network for cell line using SPIDERand message-passing from PANDA as a matrix of size (t,g).
%
% Authors: 
%               Abhijeet Sonawane, Kimberly Glass
% 
% 
% Publications:
% 
%                


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
