function [gt]=readGt(filename)
% Description:
%             Reads ground truth network in the following gmt-like format
%             tf1,gene1,gene2,
%             For some reason, the regular readtable and readmatrix
%             Functions of MATLAB fail to correctly parse this type of
%             files.
% Inputs:
%             filename: name of ground truth file
% Outputs:
%             gt: ground truth network as MATLAB table
% Example:
%             this example requires AWS CLI to download data
%             system('aws s3 cp s3://granddb/optPANDA/groundTruths/LCL_CHIP_SEQ_REMAP_ALL.txt .')
%             gtFunbind=readGt('LCL_CHIP_SEQ_REMAP_ALL.txt');
%             gtTFNamesFunbind = gtFunbind{:,1};
% Author(s):
%             Marouen Ben Guebila 01/20

    fid = fopen(filename);
    tline = fgetl(fid);
    gt=repmat({''},700,40000);
    i=0;
    lenVec=[];
    while ischar(tline)
        i=i+1;
        tline = fgetl(fid);
        if ischar(tline)
            C=strsplit(tline,',');
            gt(i,1:length(C))=C;
            lenVec=[lenVec length(C)];
        end
    end
    fclose(fid);
    gt(i:end,:)=[];
    gt(:,max(lenVec)+2:end)=[];% leave one empty space for later parsing
    gt=cell2table(gt);
    
end