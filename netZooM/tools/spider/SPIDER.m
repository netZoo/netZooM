function RegNet=SPIDER(RegNet, GeneCoReg, TFCoop, alpha);
% Description:
%             Using SPIDER to infer epigenetically-informed gene regulatory network. This function uses PANDA at the back end of the computation 
%             1. Normalizing networks
%             2. Running PANDA algorithm
%             3. Writing out SPIDER network (optional)
% Inputs:
%             alpha        : parameter that determines the level of message-passing
%		      RegNet       : motif prior of gene-TF regulatory network
%             GeneCoReg    : Optional input for SPIDER: Identify matrix of size GeneNames by GeneNames is used as input by default where gene expression information is absent
%             TFCoop       : Optional input for SPIDER: Identify matrix of size TFNames by TFNames is used as input by default if PPI network is not used as input
%             TFCoop       : PPI binding between transcription factors
% Outputs:
%             SpiderNet     : Predicted TF-gene gene complete regulatory network for cell line using SPIDER and message-passing from PANDA as a matrix of size (t,g).
% Author(s): 
%             Abhijeet Sonawane, Kimberly Glass

    [NumTFs,NumGenes]=size(RegNet);

    %% Run PANDA %%

    disp('Adjusting the degree of prior network!');
    RegNet=DegreeAdjust(RegNet);
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
    
end
