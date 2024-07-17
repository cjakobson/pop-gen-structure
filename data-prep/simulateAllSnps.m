%simulate all coding SNPs in S288c
%categorize by syn/mis and Ts/Tv
%cross-reference to secondary structures from alphafold


clear all

figureCounter=1; %figure counter

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

inputGenome=fastaread('orf_coding.fasta');

mkdir('mutationOutput')

systematicName=cell(length(inputGenome),1);


%nGenesToDo=100;
nGenesToDo=length(systematicName);


tic

for i=1:nGenesToDo
    
    tempStr=strsplit(inputGenome(i).Header,' ');
    systematicName{i}=tempStr{1};
    dnaSequence=inputGenome(i).Sequence;
    
    %add TsTv to this
    simulateSnps(systematicName{i},dnaSequence);
    
    if mod(i,100)==0
        i
    end
    
end

toc



tic

%now get secondary structure for each residue
pdbDictionary=tdfread('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/UP000002311_559292_YEAST/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.07.07-20.59.31.22.tsv');

pdbId=cell(length(pdbDictionary.Entry),1);
systematicNamePdb=pdbId;
m=1;
for i=1:length(pdbDictionary.Entry)
    
    tempStr=pdbDictionary.Entry(i,:);
    tempStr(tempStr==' ')=[];
    pdbId{m}=tempStr;
    
    tempStr=pdbDictionary.Gene_Names_0x28ordered_locus0x29(i,:);
    tempStr(tempStr==' ')=[];
    systematicNamePdb{m}=tempStr;
    m=m+1;
    
end

%pSite data
%phosphosite analysis
phosphoData=readtable('Lanz2021/embr202051121-sup-0003-datasetev2.xlsx',...
    'sheet','EV2_Lanz et al. Dataset');

for i=1:nGenesToDo
    
    if sum(ismember(systematicNamePdb,systematicName{i}))>0
        
        pdbQuery=pdbId(ismember(systematicNamePdb,systematicName{i}));
        pdbQuery=pdbQuery{1};   %in case of multiple mapping

        %include dnaSequence to refer to later wrt 1k genomes (to get
        %codons)
        dnaSequence=inputGenome(i).Sequence;
        
        %getDssp(systematicName{i},pdbQuery,dnaSequence);
        %getNeighbors(systematicName{i},pdbQuery);
        %getPsites(phosphoData,systematicName{i});
        
        if mod(i,100)==0
            
            i
            
        end
        
    end
    
end

toc



tic

%match tables and get summary statistics
%move this to another script or comment out the above once have all tables made

misSecondary=[];
synSecondary=[];

misSasa=[];
synSasa=[];

misNormSasa=[];
synNormSasa=[];

misPhi=[];
misPsi=[];

synPhi=[];
synPsi=[];

misTs=[];
synTs=[];

misNewCodon=[];
misInitialCodon=[];
synNewCodon=[];
synInitialCodon=[];

misPosInOrf=[];
misPosInRun=[];
misLengthOfRun=[];
synPosInOrf=[];
synPosInRun=[];
synLengthOfRun=[];

misNeighbors=[];
synNeighbors=[];

misPdist=[];
synPdist=[];

tic
for i=1:nGenesToDo
    
    if sum(ismember(systematicNamePdb,systematicName{i}))>0
        pdbQuery=pdbId(ismember(systematicNamePdb,systematicName{i}));
        pdbQuery=pdbQuery{1};   %in case of multiple mapping
        toGet4=['mutationOutput/' pdbQuery 'pDist.mat'];
    else
        toGet4='NA';
    end
    
    toGet1=['mutationOutput/' systematicName{i} 'mutationTable.mat'];
    toGet2=['mutationOutput/' systematicName{i} 'dsspTable.mat'];
    toGet3=['mutationOutput/' systematicName{i} 'neighborTable.mat'];
    
    if exist(toGet1)&&exist(toGet2)&&exist(toGet3)
        
        load(toGet1)
        mutationTable=outputTable;
        %clip stop codon mutants (not in DSSP)
        mutationTable(mutationTable.residueNumber==max(mutationTable.residueNumber),:)=[];
        %get rid of start codon also? not comparable to mis/syn from 1k
        
        load(toGet2)
        dsspTable=outputTable;
        
        load(toGet3)
        neighborTable=outputTable;
        
        if exist(toGet4)
            load(toGet4)
            pDistTable=outputTable;
        end
        
        vMis=mutationTable.isMis==1;
        vSyn=~vMis;
        
        tempMisTs=mutationTable.isTs(vMis);
        tempSynTs=mutationTable.isTs(vSyn);
        
        %turn ref/alt codons into doubles
        tempMisInitialCodon=codonTypes(mutationTable.initialCodon(vMis));
        tempMisNewCodon=codonTypes(mutationTable.newCodon(vMis));
        
        tempSynInitialCodon=codonTypes(mutationTable.initialCodon(vSyn));
        tempSynNewCodon=codonTypes(mutationTable.newCodon(vSyn));
        
        misResidues=mutationTable.residueNumber(vMis);
        synResidues=mutationTable.residueNumber(vSyn);
        %some proteins are different lengths (alternative start?) --
        %discard
        
        %should do explicit check of length matching (could be too
        %short, as well)

        if (max(misResidues)<=height(dsspTable))&&...
                (max(synResidues)<=height(dsspTable))

            tempMisSecondary=table2array(dsspTable(misResidues,2));
            tempMisSecondary(cellfun(@isempty,tempMisSecondary))={'NA'};
            [~,tempMisSecondary]=structureTypes(tempMisSecondary);
            tempSynSecondary=table2array(dsspTable(synResidues,2));
            tempSynSecondary(cellfun(@isempty,tempSynSecondary))={'NA'};
            [~,tempSynSecondary]=structureTypes(tempSynSecondary);
            
            tempMisPosInOrf=table2array(dsspTable(misResidues,6));
            tempMisPosInRun=table2array(dsspTable(misResidues,7));
            tempMisLengthOfRun=table2array(dsspTable(misResidues,8));

            tempSynPosInOrf=table2array(dsspTable(synResidues,6));
            tempSynPosInRun=table2array(dsspTable(synResidues,7));
            tempSynLengthOfRun=table2array(dsspTable(synResidues,8));
        

            misSecondary=[misSecondary;tempMisSecondary'];
            synSecondary=[synSecondary;tempSynSecondary'];
            
            misSasa=[misSasa;table2array(dsspTable(misResidues,3))];
            synSasa=[synSasa;table2array(dsspTable(synResidues,3))];
            
            %also normalize to sasa per-protein
            tempMean=mean(table2array(dsspTable(:,3)),'omitnan');
            
            misNormSasa=[misNormSasa;table2array(dsspTable(misResidues,3))./tempMean];
            synNormSasa=[synNormSasa;table2array(dsspTable(synResidues,3))./tempMean];
            
            misPhi=[misPhi;table2array(dsspTable(misResidues,4))];
            misPsi=[misPsi;table2array(dsspTable(misResidues,5))];
            synPhi=[synPhi;table2array(dsspTable(synResidues,4))];
            synPsi=[synPsi;table2array(dsspTable(synResidues,5))];
            
            misTs=[misTs;tempMisTs];
            synTs=[synTs;tempSynTs];
            
            misNewCodon=[misNewCodon;tempMisNewCodon];
            misInitialCodon=[misInitialCodon;tempMisInitialCodon];
            
            synNewCodon=[synNewCodon;tempSynNewCodon];
            synInitialCodon=[synInitialCodon;tempSynInitialCodon];
            
            misPosInOrf=[misPosInOrf;tempMisPosInOrf];
            misPosInRun=[misPosInRun;tempMisPosInRun];
            misLengthOfRun=[misLengthOfRun;tempMisLengthOfRun];
            
            synPosInOrf=[synPosInOrf;tempSynPosInOrf];
            synPosInRun=[synPosInRun;tempSynPosInRun];
            synLengthOfRun=[synLengthOfRun;tempSynLengthOfRun];
            
            misNeighbors=[misNeighbors;neighborTable.neighbors(misResidues)];
            synNeighbors=[synNeighbors;neighborTable.neighbors(synResidues)];
            
            if exist(toGet4)
                
                misPdist=[misPdist;pDistTable.pDist(misResidues)];
                synPdist=[synPdist;pDistTable.pDist(synResidues)];
                
            else
                
                misPdist=[misPdist;nan(size(misResidues))];
                synPdist=[synPdist;nan(size(synResidues))];
                
            end
            
        end

        if mod(i,100)==0
            
            i
            toc
            
        end
        
    end

end





%change cell arrays to double to avoid huge save files
save('simulationStructureData.mat')

% 
% save('misSecondary.mat','misSecondary')
% save('synSecondary.mat','synSecondary')
% 
% save('misSasa.mat','misSasa')
% save('synSasa.mat','synSasa')
% 
% save('misTs.mat','misTs')
% save('synTs.mat','synTs')
% 
% save('misInitialCodon.mat','misInitialCodon')
% save('misNewCodon.mat','misNewCodon')
% 
% save('synInitialCodon.mat','synInitialCodon')
% save('synNewCodon.mat','synNewCodon')
% 


