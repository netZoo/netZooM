function [Adj, TFNames, GeneNames]=BuildPriorFromRegions(motifhitfile, regfile, bedtoolspath);

% create random output file tags
rval=round(rand(1)*1000000);
rtag1=['temp', num2str(rval), '-a.txt'];
rtag2=['temp', num2str(rval), '-b.txt'];

% program parameters
% bedtoolspath=''; % set equal to '' if bedtools is already on the system path
% motifhitfile='../InputData/FilteredMotifFiles/A549_filtered_motiflocations.bed'; 
% regfile='../InputData/RegulatoryRegions/DistalRegulatoryRegions.bed';

disp('Identifying overlap between motif-hits and regulatory regions');
btag1=[bedtoolspath, 'bedtools intersect -a '];
btag=['!', btag1, regfile, ' -b ', motifhitfile, ' -wo | awk ''{print $8"\t"$4"\t"$9 > "', rtag1, '"}'''];
eval(btag);
btag2=['!sort -u ', rtag1, ' > ', rtag2];
eval(btag2);

disp('Contructing Regulatory Network');
% read in mapping information for conversion to Input Network
[TF,gene,weight]=textread(rtag2, '%s%s%f');
TFNames=unique(TF);

% identify set of Genes being mapped to
[~,~,~,GeneNames]=textread(regfile, '%s%u%u%s');
GeneNames=unique(GeneNames);

% declare initial adjacency matrix variable
Adj=zeros(length(TFNames), length(GeneNames));
[~,i]=ismember(TF, TFNames);
[~,j]=ismember(gene, GeneNames);

% sort weights in decreasing order, this will make sure highest possible weight for each edge is listed first
[~,idx]=sort(weight, 'descend');
i=i(idx); j=j(idx); weight=weight(idx);
% find the unique set of edges, the first instance is the highest weight (see above), use this to weight edge
[~,idx]=unique([i,j], 'rows');
Adj=sparse(i(idx),j(idx),weight(idx), length(TFNames), length(GeneNames));
Adj=full(Adj);

% clean-up
eval(['!rm -f ', rtag1]);
eval(['!rm -f ', rtag2]);
