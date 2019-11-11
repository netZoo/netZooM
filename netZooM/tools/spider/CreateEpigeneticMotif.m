function CreateEpigeneticMotif(epifile, motifdir, outname, bedtoolspath);


% Description:
%               Create epigenetically-filtered motif locations using bedtools and combining them into one bed file 
%
% Inputs:
%       	epifile      : path to file with open chromatin regions for given cell line 
%               motifidir    : path to file containing epigenetically informed motif information, can be created using CreateEpigeneticMotif.m

% 
% Outputs:
%               outname      : path to the directory save bedfiles containing all the epigenetically-filtered motifs in one file e.g. A549_filtered_motiflocations.bed     
%
% Authors: 
%               Abhijeet Sonawane, Kimberly Glass


% program parameters
% bedtoolspath=''; % set equal to '' if bedtools is already on the system path
% epifile='../InputData/DNaseBedFiles/A549_DnasePeaks.bed'; % file with regions of open chromatin
% motifdir='../InputData/MotifBedFiles/'; % where the original motif scan files are stored (one bed file per motif)
% outname='../InputData/FilteredMotifFiles/A549_filtered_motiflocations.bed';

% identify motif files to parse
motiffiles=dir([motifdir, '*.bed']);
NumTF=length(motiffiles);
TFNames=cell(NumTF,1);

% define code snippets
btag1=[bedtoolspath, 'bedtools intersect -a '];
btag2=[' | awk ''{print $1"\t"$2"\t"$3"\t""'];
btag3=['""\t"1.0 >> "', outname, '"}'''];

% create the new output file 
eval(['!rm -f ', outname]);
eval(['!touch ', outname]);

% go through motifs one-by-one
tic
for(cnt=1:NumTF)
	TFNames{cnt}=motiffiles(cnt).name(1:end-4);
	disp(['Now working on ', TFNames{cnt}, '......']);

	% run bedtools intersect btw the motif and dnase data
	btag=['!', btag1, motifdir, motiffiles(cnt).name, ' -b ', epifile, btag2, TFNames{cnt}, btag3];
	eval(btag);

end
% update to user
timelapse=toc;
disp(['Integrating motif and epigenetic information took ', num2str(timelapse), ' seconds.']);
