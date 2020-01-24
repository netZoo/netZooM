function [ROC,gtNet,AgNet,PR,commonTFs]=validateNetwork(AgNet,gt,gtTFNames,TFNames,GeneNames)
% Description:
%               validates a TF-gene bipartite regulation network against
%               ground truth network. 
%
% Inputs:
%               AgNet     : TF-by-gene matrix representing predicted regulation
%                           network
%               gt        : ground truth experimental network. Tf-by-gene
%                           matrix
%               gtTFNames : ground truth TF names corresponding to
%                           the order of rows of gt
%               TFNames   : TF names of predicted network corresponding to
%                           the order of rows of AgNet
%               GeneNames : ground truth gene names corresponding to
%                           the order of column of gt
%                          
% Outputs:
%               ROC   : structure containing data of ROC curve
%                       ROC.T: thresholds
%                       ROC.X: false negative rate
%                       ROC.Y: true positive rate
%                       ROC.AUC: Area under the ROC curve
%               gtNet : ground truth network with commonTFs
%               AgNet : predicted network with commonTFs
%               PR    : structure containing data of PR curve
%                       PR.T: thresholds
%                       PR.X: recall
%                       PR.Y: precision
%                       PR.AUC: Area under the PR curve
%               commonTFs : TFs that intersect both gtNet and AgNet
%
% Authors: 
%               Marouen Ben Guebila 01/2020
    

    % build ground truth network
    gtf         = gt;

    % reduce predicted and gt network 
    [commonTFs, ipredTF, igtTFs] = intersect(TFNames ,gtTFNames);

    % reduce predicted network
    AgNet          = AgNet(ipredTF,:);
    
    % build ground truth network
    gtNet = zeros(length(gtTFNames),length(GeneNames));

    for i=1:length(gtTFNames)
        [a,ib,ic]  = intersect(gtf{i,2:end}, GeneNames);
        gtNet(i,ic)= 1;
    end
    
    % /!\ reorder gtNet to have the same order as AgNet
    gtNet = gtNet(igtTFs,:); 
    
    % ROC curves
    [X,Y,T,AUC] = perfcurve(gtNet(:),AgNet(:),1);
    ROC.X  =X;
    ROC.Y  =Y;
    ROC.T  =T;
    ROC.AUC=AUC;
    
     % ROC curves
    [X,Y,T,AUC] = perfcurve(gtNet(:),AgNet(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
    PR.X  =X;
    PR.Y  =Y;
    PR.T  =T;
    PR.AUC=AUC;
end
