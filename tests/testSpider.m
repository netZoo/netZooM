function test_suite=testSpider()
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function testSpiderSimple()
	% Tell if this is Octave (Unit tests) or Matlab
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

        % Load statistics package from Octave
        if isOctave
                %we need the nan package because it has a fast implementation of corrcoeff
                %pkg load statistics
            pkg load nan;
        end

        % Set Program Parameters
        motifhitfile = 'tests/spider/output/A549_filtered_motiflocations.bed'; % file storing epigenetically informed motif information, can be created with CreateEpigeneticMotif.m
        regfile      = 'tests/spider/RegulatoryRegions_0-1kb.bed'; % file containing regulatory regions for genes, can be created with DefineRegulatoryRegions.m
        
        annofile     = 'tests/spider/refseq_hg19_05292018'; % file with gene annotations
        chrinfo      = 'tests/spider/GenomeWideRanges.bed'; % file with chromosome information
        ranges       = '';%{[-1000,+1000]};
        
        motifdir     = 'tests/spider/motifs/'; % where the original motif scan files are stored (one bed file per motif)
        epifile      = 'tests/spider/A549_DnasePeaks.bed'; % file with open chromatin regions
        bedtoolspath = './bedtools2/bin/'; 
        outtag = 'tests/output/';
        
        spider_out  = 'tests/spider/output/A549_5TF_100Genes_casenet.txt';  % optional, leave empty if file output is not required
        save_pairs  = 0;%saving in .pairs format
        save_temp   = '';  % optional, leave empty if temp data files are not needed afterward
        lib_path    = '';  % path to the folder of PANDA source code
        alpha       = 0.1;
        nTF = 5; %Number of TFs in prior

        % Add path
        addpath(genpath(fullfile(pwd,'tests')));
        
        if 0
            % Call SPIDER
            SpiderNet = spider_run(lib_path, bedtoolspath, alpha, motifhitfile,  annofile,...
                chrinfo, ranges, regfile, outtag,motifdir, epifile,save_temp,save_pairs,spider_out,nTF );
        end
        
        % Now try the step-by-step approach
        CreateEpigeneticMotif(epifile, motifdir, motifhitfile, bedtoolspath, nTF);
        % Build SPIDER prior
        [PriorNet, TFNames, GeneNames]=BuildSPIDERprior(motifhitfile, regfile, bedtoolspath);
        %temporary reduction in number of genes for Actions
        numGenes = 100; 
        GeneNames= GeneNames(1:numGenes);
        PriorNet = PriorNet(:,1:numGenes);
        SpiderNet= SPIDER(PriorNet, eye(length(GeneNames)), eye(length(TFNames)), alpha);
        
        % Load the expected result
        ExpSpiderNet = textread('tests/spider/output/A549_5TF_100Genes_testnet.txt');%different behavior with Octave and Matlab
        % /!\ ExpAgNet is a row-major matrix, while reshape transforms in column-major format, thus the transpose
        ExpSpiderNet = reshape(ExpSpiderNet,[size(SpiderNet,1), size(SpiderNet,2)]);
        
        % Compare the outputs
        tolMat=1e-6;
        deltaMat=max(max(abs(SpiderNet-ExpSpiderNet)))
        assertTrue(deltaMat < tolMat);
end
