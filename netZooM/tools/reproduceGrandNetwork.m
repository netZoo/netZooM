netZoo     = 'netZooM';
netZooVer  = '0.5.1';
netZooTool = 'OTTER';
ppiLink    = 'https://granddb.s3.amazonaws.com/cancer/cervix_cancer/cancer_cervix_otter_ppi.csv';
motifLink  = 'https://granddb.s3.amazonaws.com/cancer/cervix_cancer/cancer_cervix_otter_motif.csv';
expLink    = 'https://granddb.s3.amazonaws.com/cancer/cervix_cancer/cancer_cervix_expression_tcga.csv';
networkLink= 'https://granddb.s3.amazonaws.com/cancer/cervix_cancer/cancer_cervix_otter_network.csv';
netRes=reproduceGrandNetworkFunction(ppiLink,motifLink,expLink,netZoo,netZooTool,...
    netZooVer,networkLink);

function [netRes,repPres]=reproduceGrandNetworkFunction(ppiLink,motifLink,expLink,...
    netZoo,netZooTool,netZooVer,networkLink)
% Description:
%               reproduceGrandNetworkFunction fetches input data from GRAND (grand.networmedicine.org) 
%               to recompute a regulatory network and compare it to the
%               network from GRAND. Thsi function checks that GRAND
%               networks are reproducible using netZoo (netzoo.github.io).
%
% Inputs:
%               ppiLink    : url to file containing TF-TF interaction graph
%               motifLink  : url to file containing the prior regulatory network
%               expLink    : url to file containing gene expression as a
%                            matrix of size (g,s) with s number of gene expression
%                            samples and g the number of genes
%               netZoo     : which implementation of netZoo e.g., 'netZooM'
%               netZooTool : which zoo animal e.g., 'PUMA'
%               netZooVer  : which netZoo version e.g., 0.5.1
%               networkLink: url to file containing GRAND regulatory
%                            network to comapre to
% 
% Outputs:
%               netRes     : Predicted bipartite gene complete regulatory network 
%                            using netZooTool as a matrix of size (t,g) with t the number of regulators.
%               repPres    : largest absolute difference between recomputed network
%                             and GRAND network
%
% Dependencies:
%               AWSCLI    : https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html
%               GIT       : https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
%
% Authors: 
%               Marouen Ben Guebila 6/2020

    % Parse args
    if isequal(netZooTool,'PANDA-LIONESS')
        netZooTool = 'LIONESS';
    end
    % Check dependencies
    % 1. AWSCLI
    try
       system('aws --version') 
    catch
        error('Please install awscli (https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) and try again.')
    end
    % 2. GIT
    try
       system('git --version') 
    catch
        error('Please install git (https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and try again.')
    end

    % Fetch netZoo
    pathToZoo= fileparts([which([netZooTool '.m']) which([netZooTool '_run.m'])]);
    if isempty(pathToZoo)
        system(['git clone https://github.com/netZoo/' netZoo '.git']);
        % Add to path
        addpath(genpath(netZoo))
        % Go to the root of netZoo
        eval(['cd ' netZoo]);
    else
        % Go to the root of netZoo
        cd(pathToZoo);cd ../..
    end
    % Switch to specific version
    system(['git checkout tags/' netZooVer]);

    % Fetch priors
    % change prefix to ssh
    motifLinkssh   =['s3://granddb/' motifLink(34:end)];
    expLinkssh     =['s3://granddb/' expLink(34:end)];
    ppiLinkssh     =['s3://granddb/' ppiLink(34:end)];
    networkLinkssh =['s3://granddb/' networkLink(34:end)];
    % Motif
    system(['aws s3 cp ' motifLinkssh ' .']);
    % Exp
    system(['aws s3 cp ' expLinkssh ' .']);
    % PPI
    if ~isequal(ppiLinkssh,'-')
        system(['aws s3 cp ' ppiLinkssh ' .'])
    end
    % GRAND network for comparison
    system(['aws s3 cp ' networkLinkssh ' .'])

    % Run tool
    % parse file name
    % Network
    [~,fileNet,extNet]=fileparts(networkLink);
    net_file  =[fileNet extNet];
    % Expression
    [~,fileExp,extExp]=fileparts(expLink);
    exp_file  =[fileExp extExp];
    % Motif
    [~,fileMotif,extMotif]=fileparts(motifLink);
    motif_file=[fileMotif extMotif];
    % PPI
    [~,filePPI,extPPI]=fileparts(ppiLink);
    ppi_file  =[filePPI extPPI];
    lib_path='';
    switch netZooTool
        case 'PANDA'
            netRes=panda_run(lib_path, exp_file, motif_file, ppi_file);
        case 'LIONESS'    
            lib_path='';
            netRes=lioness_run(lib_path, exp_file, motif_file, ppi_file);
        case 'PUMA'   
            outtag   = '';
            alpha    = 0.1;
            mir_file = '';
            netRes=RunPUMA(outtag,alpha,motif_file,exp_file,ppi_file,mir_file);
        case 'OTTER'
            PPI  = readtable(ppi_file,'ReadVariableNames',1,'ReadRowNames',1,...
                'PreserveVariableNames',1);
            Expr = readtable(exp_file,'ReadVariableNames',1,'ReadRowNames',1,...
                'PreserveVariableNames',1);
            Motif= readtable(motif_file,'ReadVariableNames',1,'ReadRowNames',1,...
                'PreserveVariableNames',1);
            C = corrcoef(Expr{:,:}');
            netRes=otter(Motif{:,:},PPI{:,:},C);
    end

    % Compare to result network
    netGrand  = readtable(net_file,'ReadVariableNames',1,'ReadRowNames',1);
    % compute precision difference to check reproducibility
    repPres = max(max( abs(netGrand{:,:} - netRes) ));

end