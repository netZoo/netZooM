outtag='PUMA';

% Set Program Parameters
alpha=0.1;
motif_file='ToyData/ToyMotifData.txt';
exp_file='ToyData/ToyExpressionData.txt';
ppi_file='ToyData/ToyPPI.txt';
% ppi_file=''; % PUMA/PANDA can be run without PPI data
% mir_file=''; % this will run PANDA
mir_file='ToyData/ToyMiRList.txt'; % run PUMA

% select samples
% if you want to select for example samples 201-300, set SelectSize=100 and Offset=200
SelectSize=20; % nr of samples to run
Offset=0; % offset


%% Read in Data %%
disp('Reading in data!')

% Expression Data
fid=fopen(exp_file, 'r');
headings=fgetl(fid);
NumConditions=length(regexp(headings, '\t'));
frewind(fid);
Exp=textscan(fid, ['%s', repmat('%f', 1, NumConditions)], 'delimiter', '\t', 'CommentStyle', '#');
fclose(fid);
GeneNames=Exp{1};
NumGenes=length(GeneNames);
Exp=cat(2, Exp{2:end});
% Exp=quantilenorm(Exp); % quantile normalize data

% Prior Regulatory Network and miR information (PUMA)
[TF, gene, weight]=textread(motif_file, '%s%s%f');
if(~isempty(mir_file)) % check miRs in mir_file
	[miR]=textread(mir_file, '%s');
end

TFNames=unique(TF);
NumTFs=length(TFNames);
[~,i]=ismember(TF, TFNames); % convert regulators to integers
[~,j]=ismember(gene, GeneNames); % convert genes to integers
RegNet=zeros(NumTFs, NumGenes);
RegNet(sub2ind([NumTFs, NumGenes], i, j))=weight; % put weights from motif_file in matrix

% if mir_file exists (PUMA), check indices of miRs in the PPI data
if(~isempty(mir_file))
	[~,k]=ismember(miR, TFNames);
        m=(1:NumTFs)';
        [s1,s2] = ndgrid(k,m);
        [t1,t2] = ndgrid(m,k);
end

% PPI Data
TFCoop=eye(NumTFs); % identity matrix of PPI file, or ppi_file
if(~isempty(ppi_file))
	[TF1,TF2,weight]=textread(ppi_file, '%s%s%f');
	[~,i]=ismember(TF1, TFNames);
	[~,j]=ismember(TF2, TFNames);
	TFCoop(sub2ind([NumTFs, NumTFs], i, j))=weight;
	TFCoop(sub2ind([NumTFs, NumTFs], j, i))=weight;
	% set all edges between regulators present in mirlist and other regulars 0
	if(~isempty(mir_file))
		TFCoop(sub2ind([NumTFs, NumTFs],s1,s2))=zeros(size(k,1), NumTFs);
		TFCoop(sub2ind([NumTFs, NumTFs],t1,t2))=zeros(NumTFs, size(k,1));
	end
	% set diagonal of TFCoop to 1 (self-interactions)
	TFCoop(1:(size(TFCoop,1)+1):end) = 1; % used in both PANDA and PUMA
end

NumConditions=size(Exp,2);

GeneCoReg=corr(Exp', 'type', 'pearson', 'rows', 'pairwise');
%NumUsed=double(~isnan(Exp))*double(~isnan(Exp)'); % optional: determine nr of samples with expression values (not NaN)
%GeneCoReg=GeneCoReg.*(NumUsed/NumConditions); % optional: normalize the correlation based on NumUsed
GeneCoReg(1:NumGenes+1:NumGenes^2)=1; % set diagonal to 1
GeneCoReg(isnan(GeneCoReg))=0; % change NAs to 0

if(isempty(mir_file)) % run PANDA
	AgNet=PANDA(RegNet, GeneCoReg, TFCoop, alpha);
	AgNet=AgNet(:);
end
if(~isempty(mir_file)) % run PUMA
	AgNet=PUMA(RegNet, GeneCoReg, TFCoop, alpha, s1, s2, t1, t2); % s1, s2, t1, t2 are indices of miR interactions in TFCoop
	AgNet=AgNet(:);
end

% this code will run LIONESS to generate single-sample networks
myNumConditions=(Offset+1):(Offset+SelectSize);
PredNet=zeros(NumTFs*NumGenes, SelectSize);
for(condcnt=myNumConditions)
	idx=[1:(condcnt-1),(condcnt+1):NumConditions];
	GeneCoReg=corr(Exp(:,idx)', 'type', 'pearson', 'rows', 'pairwise');
	%NumUsed=double(~isnan(Exp(:,idx)))*double(~isnan(Exp(:,idx))');  % optional: determine nr of samples with expression values (not NaN)
	%GeneCoReg=GeneCoReg.*(NumUsed/(NumConditions-1)); % optional: normalize the correlation based on NumUsed
	GeneCoReg(1:NumGenes+1:NumGenes^2)=1; % set diagonal to 1
	GeneCoReg(isnan(GeneCoReg))=0; % change NAs to 0
	if(isempty(mir_file)) % run PANDA
		LocNet=PANDA(RegNet, GeneCoReg, TFCoop, alpha);
	end
	if(~isempty(mir_file)) % run PUMA
		LocNet=PUMA(RegNet, GeneCoReg, TFCoop, alpha, s1, s2, t1, t2);
	end
	idPredNet=condcnt-Offset;
	PredNet(:,idPredNet)=NumConditions*(AgNet-LocNet(:))+LocNet(:);
end


TF=repmat(TFNames, 1, length(GeneNames));
gene=repmat(GeneNames', length(TFNames), 1);
TF=TF(:);
gene=gene(:);
RegNet=RegNet(:);

fid=fopen([outtag, '_FinalNetwork.pairs'], 'wt');
for(cnt=1:length(TF))
	fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, RegNet(cnt), AgNet(cnt));
end
fclose(fid);


% to save the LIONESS networks in a .mat file
save([outtag, '_LIONESSNetworks.mat'], '-v7.3', 'PredNet', 'TF', 'gene', 'AgNet', 'RegNet');

% optional: print single-sample edge weights in a .txt file
%dlmwrite([outtag,'_LIONESSNetworks.txt'], PredNet, '\t');