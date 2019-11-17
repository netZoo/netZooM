
% Set Program Parameters

% Set Program Parameters
motifhitfile = 'tests/spider_data/output/A549_filtered_motiflocations.bed'; % file storing epigenetically informed motif information, can be created with CreateEpigeneticMotif.m
regfile     = 'tests/spider_data/RegulatoryRegions_0-1kb.bed'; % file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m

annofile    = 'tests/spider_data/refseq_hg19_05292018'; % file with gene annotations
chrinfo     = 'tests/spider_data/GenomeWideRanges.bed'; % file with chromosome information
ranges      = ''; %Sample{[-1000,+500]}
motifdir    = 'tests/spider_data/motifs/'; % where the original motif scan files are stored (one bed file per motif)
epifile     = 'tests/spider_data/A549_DnasePeaks.bed'; % file with open chromatin regions
bedtoolspath = ''  %to be specified by Marouen
outtag = 'tests/output/';
spider_out  = 'tests/spider_data/output/A549_5TF_testnet.txt';  % optional, leave empty if file output is not required
save_pairs = 0;%saving in .pairs format
save_temp  = '';  % optional, leave empty if temp data files are not needed afterward
lib_path   = '/udd/spaso/netZooM';  % path to the folder of PANDA source code
alpha      = 0.1;
nTF = 5; %Number of TFs in prior


