function RegNet=PUMA(RegNet, GeneCoReg, TFCoop, alpha, s1, s2, t1, t2);
% Description:
%               PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators 
%               of mRNA expression levels. 
%
% Inputs:
%               RegNet   : motif prior of gene-TF regulatory network
%               GeneCoReg: gene-gene co-regulatory network
%               TFCoop   : PPI binding between transcription factors
%               alpha    : learning rate
%               s1, s2, t1, t2 are indices of miR interactions in TFCoop, the TF-TF PPI network.
%
% Outputs:
%               RegNet   : inferred gene-TF regulatory network
%
% Authors:
%               Marieke Kuijjer
%
% Publications:
%               https://www.ncbi.nlm.nih.gov/pubmed/28506242

[NumTFs,NumGenes]=size(RegNet);

%% Run PANDA %%

disp('Normalizing Networks!');
RegNet=NormalizeNetwork(RegNet);
GeneCoReg=NormalizeNetwork(GeneCoReg);
TFCoop=NormalizeNetwork(TFCoop);
TFCoopInit=TFCoop;% PUMA, keep a backup of the initial TFCoop
% you should keep those values for the ismember(miR, TFNames)

tic;
disp('Learning Network!')
step=0;
hamming=1;
while(hamming>0.001)
	Responsibility=Tfunction(TFCoop, RegNet);
	Availability=Tfunction(RegNet, GeneCoReg);
	hamming=sum(abs(RegNet(:)-0.5*(Responsibility(:)+Availability(:))))/(NumTFs*NumGenes);
	RegNet=(1-alpha)*RegNet+alpha*0.5*(Responsibility+Availability);

	PPI=Tfunction(RegNet, RegNet');
	PPI=UpdateDiagonal(PPI, NumTFs, alpha, step);
	TFCoop=(1-alpha)*TFCoop+alpha*PPI;
	TFCoopDiag=diag(TFCoop);
	% PUMA
	TFCoop(sub2ind([NumTFs, NumTFs],s1,s2))=TFCoopInit(sub2ind([NumTFs, NumTFs],s1,s2)); % PUMA
	TFCoop(sub2ind([NumTFs, NumTFs],t1,t2))=TFCoopInit(sub2ind([NumTFs, NumTFs],t1,t2)); % PUMA
        TFCoop(1:(size(TFCoop,1)+1):end) =TFCoopDiag; % PUMA

	CoReg2=Tfunction(RegNet', RegNet);
	CoReg2=UpdateDiagonal(CoReg2, NumGenes, alpha, step);
	GeneCoReg=(1-alpha)*GeneCoReg+alpha*CoReg2;

	disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
	step=step+1;
end
runtime=toc;
NummiRs=size(unique(s1),1);
fprintf('Running PUMA on %d Genes and %d regulators (%d TFs and %d miRs) took %f seconds!\n', NumGenes, NumTFs, NumTFs-NummiRs, NummiRs, runtime);
