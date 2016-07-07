function RegNet=PANDA(RegNet, GeneCoReg, TFCoop, alpha);

[NumTFs,NumGenes]=size(RegNet);

%% Run PANDA %%

disp('Normalizing Networks!');
RegNet=NormalizeNetwork(RegNet);
GeneCoReg=NormalizeNetwork(GeneCoReg);
TFCoop=NormalizeNetwork(TFCoop);

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

	CoReg2=Tfunction(RegNet', RegNet);
	CoReg2=UpdateDiagonal(CoReg2, NumGenes, alpha, step);
	GeneCoReg=(1-alpha)*GeneCoReg+alpha*CoReg2;

	disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
	step=step+1;
end
runtime=toc;
fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
