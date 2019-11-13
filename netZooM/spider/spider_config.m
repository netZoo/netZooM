% Set Program Parameters
motifhitfile = 'X/udd/spaso/netZooM/tests/output/A549_filtered_motiflocations.bed'; % file storing epigenetically informed motif information, can be created with CreateEpigeneticMotif.m
regfile     = '/udd/spaso/netZooM/tests/spider_data/RegulatoryRegions_0-1kb.bed'; % file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m

annofile    = '/udd/spaso/netZooM/tests/spider_data/refseq_hg19_05292018'; % file with gene annotations
chrinfo     = '/udd/spaso/netZooM/tests/spider_data/GenomeWideRanges.bed'; % file with chromosome information
ranges      = {[-1000,+1000]};

motifdir    = '/udd/spaso/netZooM/tests/spider_data/motifs/'; % where the original motif scan files are stored (one bed file per motif)
epifile     = '/udd/spaso/netZooM/tests/spider_data/A549_DnasePeaks.bed'; % file with open chromatin regions
bedtoolspath = ''  %to be specified by Marouen

spider_out  = '/udd/spaso/netZooM/tests/spider.mat';  % optional, leave empty if file output is not required
save_temp  = '/udd/spaso/netZooM/tests/output/';  % optional, leave empty if temp data files are not needed afterward
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;
