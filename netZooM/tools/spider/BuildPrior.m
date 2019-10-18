function [Adj, TFNames, GeneNames]=BuildPrior(epifile, regfile, motifdir, bedtoolspath);

% create random output file tags
rval=round(rand(1)*1000000);
rtag1=['temp', num2str(rval), '-a.txt'];
rtag2=['temp', num2str(rval), '-b.txt'];

% program parameters
% bedtoolspath=''; % set equal to '' if bedtools is already on the system path
% epifile='/udd/spaso/SPIDER_Sep-17/FOLDER1_RAWDATA/DNase_data/A549_DnasePeaks.bed'; % file with regions of open chromatin
% motifdir='InputData/MotifBedFiles/'; % where the original motif scan files are stored (one bed file per motif)
% regfile='DistalRegulatoryRegions.bed';

% identify set of Genes being mapped to
[~,~,~,GeneNames]=textread(regfile, '%s%u%u%s');
GeneNames=unique(GeneNames);

% identify motif files to parse
motiffiles=dir([motifdir, '*.bed']);
NumTFs=length(motiffiles);
TFNames=cell(NumTFs,1);

% declare initial adjacency matrix variable
Adj=zeros(NumTFs, length(GeneNames));

% define code snippets
btag1=[bedtoolspath, 'bedtools intersect -a '];
btag2=[' | awk ''{print $1"\t"$2"\t"$3 > "', rtag1, '"}'''];
btag3=[' | awk ''{print $4 > "', rtag2, '"}'''];

% make sure the files to be written aren't there (in the random chance the initial outputs are empty and the files can't be written)
eval(['!rm -f ', rtag1]);
eval(['!rm -f ', rtag2]);

% go through motifs one-by-one
tic
for(cnt=1:NumTFs)
	TFNames{cnt}=motiffiles(cnt).name(1:end-4);
	disp(['Now working on ', TFNames{cnt}, '......']);

	% run bedtools intersect btw the motif and dnase data
	btag=['!', btag1, motifdir, motiffiles(cnt).name, ' -b ', epifile, btag2];
	eval(btag);

	% run bedtools intersect btw the regulatory regions and filtered motif data
	btag=['!', btag1, regfile, ' -b ', rtag1, btag3];
	eval(btag);

	if(exist(rtag2))
		% read in the file with intersected motif positions
		locG=unique(textread(rtag2, '%s'));
		Adj(cnt,:)=Adj(cnt,:)+ismember(GeneNames, locG)';
	end

	% clean-up
	eval(['!rm -f ', rtag1]);
	eval(['!rm -f ', rtag2]);
end
timelapse=toc;
disp(['Integrating motif and epigenetic information took ', num2str(timelapse), ' seconds.']);
