function RegNet=RunPUMA(outtag,alpha,motif_file,exp_file,ppi_file,mir_file,...
                        computing)
% Description:
%             PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators of mRNA expression levels.
%             The script must be run with a regulatory prior (variable `motif_file` in the RunPUMA script) and expression data (variable `exp_file`). 
%             Protein-protein interaction data (variable `ppi_file`) and a list of microRNAs (variable `mir_file`) are optional parameters. 
%             If no `mir_file` is listed, the script will run the message passing algorithm described in PANDA to estimate regulatory edges, assuming that all regulators are transcription factors. 
%             If a `mir_file` is listed, the script will run PUMA to estimate regulatory edges. In particular, regulators listed in that file will be treated as regulators that cannot form complexes with other regulators, 
%             such as microRNAs, while regulators not listed will be treated as regulators that can form complexes, such as transcription factors.
%             Some important notes on this tool:
%             - The target genes in the regulatory prior (`motif_file`) should match the target genes in the expression data (`exp_file`).
%             - The regulators in `mir_list` should match symbols in the regulatory prior (`motif_file`).
%             - Self-interactions in the protein-protein interaction data will be set to 1. The algorithm converges around these edges. Protein-protein interaction "prior" edges therefore should have weights between 0-1.
%             - In the protein-protein interaction prior, edges between microRNAs from the `mir_file` (see `RunPUMA.m` script) and other regulators can have values, 
%             but these will automatically be set to 0, as the algorithm assumes that any regulator listed in `mir_file` will not be able to form edges with other regulators.
%             Example files can be found in the folder tests/test_data/PUMA_ToyData. 
%             Some useful scripts can be found in the folder tools:
%             - bash scripts `RunPUMA.sh` and `RunPUMALIONESS.sh` can be used to remotely run `RunPUMA.m` and `RunPUMALIONESS.m`, respectively.
%             - `getCompleteEdgelist.R` is an R script that can convert an unweighted regulatory prior in a complete edgelist. This can be useful when preparing the regulatory prior and expression data.
% Inputs:
%             exp_file  : path to file containing gene expression as a matrix of size (g,g)
%             motif_file: path to file containing the prior TF-gene regulatory network based on TF motifs as a matrix of size (t,g)
%             ppi_file  : path to file containing TF-TF interaction graph as a matrix of size (t,t)
%             outtag    : path to save output PUMA network in .pairs format.
%             mir_file  : path to file containing microRNA file
%             alpha     : learning parameter for the PUMA algorithm
%             computing : 'cpu'(default)
%                         'gpu' uses GPU to compute distances
% Outputs:
%             RegNet    : Predicted TF-gene gene complete regulatory network using PANDA as a matrix of size (t,g).
% Author(s):
%             Marieke Kuijjer
% Publications:
%             https://doi.org/10.1093/bioinformatics/btaa571
%             https://www.ncbi.nlm.nih.gov/pubmed/28506242

    %% Read in Data %%
    if nargin<7
        computing='cpu';
    end

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
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        GeneCoReg=corrcoef(Exp', 'Mode', 'Pearson', 'rows', 'pairwise'); 
    else
        GeneCoReg = corr(Exp', 'type', 'pearson', 'rows', 'pairwise');
    end
    %NumUsed=double(~isnan(Exp))*double(~isnan(Exp)'); % optional: determine nr of samples with expression values (not NaN)
    %GeneCoReg=GeneCoReg.*(NumUsed/NumConditions); % optional: normalize the correlation based on NumUsed
    GeneCoReg(1:NumGenes+1:NumGenes^2)=1; % set diagonal to 1
    GeneCoReg(isnan(GeneCoReg))=0; % change NaNs to 0

    if(isempty(mir_file)) % run PANDA
        respWeight=0.5;
        similarityMetric='Tfunction';
        AgNet=PANDA(RegNet, GeneCoReg, TFCoop, alpha, respWeight,...
                    similarityMetric, computing);
        AgNet=AgNet(:);
    end
    if(~isempty(mir_file)) % run PUMA
        AgNet=PUMA(RegNet, GeneCoReg, TFCoop, alpha, s1, s2, t1, t2,...
            computing); % s1, s2, t1, t2 are indices of miR interactions in TFCoop
        AgNet=AgNet(:);
    end

    TF=repmat(TFNames, 1, length(GeneNames));
    gene=repmat(GeneNames', length(TFNames), 1);
    TF=TF(:);
    gene=gene(:);
    RegNet=RegNet(:);

    % print the network
    fid=fopen([outtag, '_FinalNetwork.pairs'], 'wt');
    for(cnt=1:length(TF))
        fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, RegNet(cnt), AgNet(cnt));
    end
    fclose(fid);
    
end
