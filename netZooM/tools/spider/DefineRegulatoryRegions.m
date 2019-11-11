function AnnoData=DefineRegulatoryRegions(annofile, ranges, outname, chrinfo);



% Description:
%               Using gene annotation and gene coordinates, define the window of user-defined ranges around TSS
%
% 	Inputs:
%               bedtoolspath : path of the bedtools (can be installed from : "https://bedtools.readthedocs.io/en/latest/content/installation.html")
%               annofile     : file with gene annotations e.g., './ReferenceData/refseq_hg19_05292018'
%               chrinfo      : file with chromosome information e.g.,  '/InputData/ReferenceData/GenomeWideRanges.bed'
%               ranges       : user-defined input for ranges around TSS for constructing proximal or distal SPIDER networks
%   	Output:
%               AnnoData     : path to save bedfile of regulatory regions used asinput to build SPIDER prior network. 

% Authors: 
%               Abhijeet Sonawane, Kimberly Glass
% 
% 
% Publications:
% 




% currently built in, could be modified
% annofile='./ReferenceData/refseq_hg19_05292018'; % file with gene annotations
% ranges={[-100000,-25000], [25000, 100000]};
% outname='DistalRegulatoryRegions.bed';
% chrinfo='../SPIDER/DNaseBedFiles/GenomeWideRanges.bed';

% read in annotation info
disp('Parsing Annotation Information')
fid=fopen(annofile, 'r');
RefGene=textscan(fid, '%u%s%s%s%f%f%f%f%f%s%s%s%s%s%s%s', 'delimiter', '\t', 'headerlines', 1);
fclose(fid);

% only keep gene annotations in the "canonical" chromosomal locations (these were the only chromosomes scanned anyway)
[allchr, chrS, chrE]=textread(chrinfo, '%s%f%f', 'delimiter', '\t');
filter=ismember(RefGene{3}, allchr);
for(cnt=1:length(RefGene))
	RefGene{cnt}=RefGene{cnt}(filter);
end

% parse annotation information for fast computing
chridx=cell(length(allchr),1);
chrTSS=cell(length(allchr),1);
chrStrand=cell(length(allchr),1);
chrStrandBool=cell(length(allchr),1);
for(cnt=1:length(allchr))
        chridx{cnt}=find(strcmp(RefGene{3}, allchr{cnt}));
        locStrand=RefGene{4}(chridx{cnt});
        locEnd=RefGene{6}(chridx{cnt});
        chrTSS{cnt}=RefGene{5}(chridx{cnt});
        chrTSS{cnt}(strcmp(locStrand, '-'))=locEnd(strcmp(locStrand, '-'));
        chrTSS{cnt}=double(chrTSS{cnt});
        chrStrand{cnt}=locStrand;

	[uStart,sidx,sloc]=unique(chrTSS{cnt});
	locGene=RefGene{13}(chridx{cnt});
	midx=zeros(100,1);
	mcnt=0;
	for(scnt=1:length(uStart))
		% if the same TSS can be attributed to multiple genes
		if(length(unique(locGene(sloc==scnt)))>1)
			lidx=find(sloc==scnt);
			[mG,fidx]=unique(locGene(lidx));
			fidx=sort(fidx, 'ascend');
			if(sidx(scnt)~=lidx(fidx(1)))
				disp('assumption violation. please type "dbquit" and notify the author');
				keyboard;
			end
			% add those genes to an vector containing the index of their locations
			for(c=2:length(mG))
				mcnt=mcnt+1;
				midx(mcnt)=lidx(fidx(c));
			end;
		end
	end
	sidx=[sidx; midx(1:mcnt)];
	chridx{cnt}=chridx{cnt}(sidx);
	chrTSS{cnt}=chrTSS{cnt}(sidx);
	chrStrand{cnt}=chrStrand{cnt}(sidx);
        chrStrandBool{cnt}=ones(length(chrStrand{cnt}),1);
	chrStrandBool{cnt}(strcmp(chrStrand{cnt}, '-'))=-1*chrStrandBool{cnt}(strcmp(chrStrand{cnt}, '-'));
end
RefGene=RefGene{13};

% collapse data into struct
disp('Collapsing Data');
NumR=length(ranges);
AnnoData=struct('chr', cell(NumR,1), 'Start', cell(NumR,1), 'End', cell(NumR,1), 'GeneName', cell(NumR,1));

for(rcnt=1:NumR)
	AnnoData(rcnt).chr=cell(length(allchr),1);
	AnnoData(rcnt).Start=cell(length(allchr),1);
	AnnoData(rcnt).End=cell(length(allchr),1);
	AnnoData(rcnt).Gene=cell(length(allchr),1);
	for(chrcnt=1:length(allchr))
		NumTSS=length(chrTSS{chrcnt});
		AnnoData(rcnt).chr{chrcnt}=repmat(allchr(chrcnt), NumTSS,1);
		AnnoData(rcnt).Start{chrcnt}=chrTSS{chrcnt};
		AnnoData(rcnt).End{chrcnt}=chrTSS{chrcnt};
		AnnoData(rcnt).Gene{chrcnt}=RefGene(chridx{chrcnt});

		sfilter=chrStrandBool{chrcnt};
		% for genes on the positive strand, range is treated "normally"
		AnnoData(rcnt).Start{chrcnt}(sfilter==1)=AnnoData(rcnt).Start{chrcnt}(sfilter==1)+ranges{rcnt}(1);
		AnnoData(rcnt).End{chrcnt}(sfilter==1)=AnnoData(rcnt).End{chrcnt}(sfilter==1)+ranges{rcnt}(2);
		% for genes on the negative strand, range is treated "opposite"
		AnnoData(rcnt).Start{chrcnt}(sfilter==-1)=AnnoData(rcnt).Start{chrcnt}(sfilter==-1)-ranges{rcnt}(2);
		AnnoData(rcnt).End{chrcnt}(sfilter==-1)=AnnoData(rcnt).End{chrcnt}(sfilter==-1)-ranges{rcnt}(1);

		% truncate ranges that fall off the edge of the chromosome
		f=AnnoData(rcnt).Start{chrcnt}<chrS(chrcnt);
		AnnoData(rcnt).Start{chrcnt}(f)=chrS(chrcnt);
		f=AnnoData(rcnt).End{chrcnt}>chrE(chrcnt);
		AnnoData(rcnt).End{chrcnt}(f)=chrE(chrcnt);

		% remove any regions that are completely off the chromosome
		f=~(AnnoData(rcnt).End{chrcnt}<chrS(chrcnt) | AnnoData(rcnt).Start{chrcnt}>chrE(chrcnt));
		AnnoData(rcnt).chr{chrcnt}=AnnoData(rcnt).chr{chrcnt}(f);
		AnnoData(rcnt).Start{chrcnt}=AnnoData(rcnt).Start{chrcnt}(f);
		AnnoData(rcnt).End{chrcnt}=AnnoData(rcnt).End{chrcnt}(f);
		AnnoData(rcnt).Gene{chrcnt}=AnnoData(rcnt).Gene{chrcnt}(f);
	end
	AnnoData(rcnt).chr=cat(1, AnnoData(rcnt).chr{:});
	AnnoData(rcnt).Start=cat(1, AnnoData(rcnt).Start{:});
	AnnoData(rcnt).End=cat(1, AnnoData(rcnt).End{:});
	AnnoData(rcnt).Gene=cat(1, AnnoData(rcnt).Gene{:});
end
AnnoData={cat(1, AnnoData(:).chr), cat(1, AnnoData(:).Start), cat(1, AnnoData(:).End), cat(1, AnnoData(:).Gene)};

% write data to file
fid=fopen(outname, 'wt');
for(cnt=1:length(AnnoData{1}))
	fprintf(fid, '%s\t%u\t%u\t%s\n', AnnoData{1}{cnt}, AnnoData{2}(cnt), AnnoData{3}(cnt), AnnoData{4}{cnt});
end
fclose(fid);
