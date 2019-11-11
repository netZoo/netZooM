% Set Program Parameters
motifhitfile = '/SPIDER_code/InputData/EpiMotifFiles/GM12878_filtered_motiflocations.bed'; % file storing epigenetically informed motif information, can be created with CreateEpigeneticMotif.m
regfile     = '/SPIDER_code/InputData/RegulatoryRegions/RegulatoryRegions_0-1kb.bed'; % file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m

annofile    = 'InputData/ReferenceData/refseq_hg19_05292018'; % file with gene annotations
chrinfo     = 'InputData/ReferenceData/GenomeWideRanges.bed'; % file with chromosome information
ranges      = {[-1000,+1000]};

motifdir    = '/SPIDER_code/InputData/MotifBedFiles/'; % where the original motif scan files are stored (one bed file per motif)
epifile     = '/SPIDER_code/InputData/DNaseBedFiles/GM12878_DnasePeaks.bed'; % file with open chromatin regions
bedtoolspath = ''

spider_out  = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/panda.mat';  % optional, leave empty if file output is not required
save_temp  = '/ifs/labs/cccb/projects/share/GTEx/LIONESS/mat/';  % optional, leave empty if temp data files are not needed afterward
lib_path   = 'lib';  % path to the folder of PANDA source code
alpha      = 0.1;
