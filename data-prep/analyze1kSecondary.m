
%assess secondary structure etc for 1k genomes variants

clear

tic

figureCounter=1; %figure counter

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

load('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/1002data/1002dataAnn.mat')

tempIdx=cellfun(@isempty,gene);

%only need gene, position, variant type, and AF for now
gene(tempIdx)=[];
proteinEncoded(tempIdx)=[];
type(tempIdx)=[];
af(tempIdx)=[];

%minor allele
af(af>0.5)=1-af(af>0.5);
af(af==0)=nan;


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


structureLabels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
    'bend','turn','unstr.'};
%structureLabels={'helix','sheet','turn','unstr'};

% 
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,3,1)
% hold on
% 
% histogram(log10(af(ismember(type,'missense_variant'))),-4:0.4:0)
% hold on
% histogram(log10(af(ismember(type,'synonymous_variant'))),-4:0.4:0)
% 
% ylabel('frequency')
% xlabel('minor allele frequency')
% legend({'mis','syn'})
% title('in 1K genomes')

essentialList=readtable('Giaever_G_et_al_2002_phenotypes_filtered_by_inv.txt');
essentialList=table2array(essentialList(:,2));
isEssential=zeros(length(gene),1);
%check essentiality
isEssential=ismember(gene,essentialList);



%mean TPMs to subset by expression level
rnaSeqData=readtable('tpmSummary.csv');
rnaSeqMat=table2array(rnaSeqData(:,2:end));

sampleInfo=readtable('sampleInfo.txt');

geneNames=rnaSeqData.ORF;

strains={'RM11','YJM975'};
conditions={'untreated','radicicol'};

meanData=[];
m=1;
for i=1:length(strains)
    for j=1:length(conditions)
        
        tempIdx=logical(ismember(sampleInfo.strain,strains{i}).*...
            ismember(sampleInfo.condition,conditions{j}));
        
        meanData{m}=mean(rnaSeqMat(:,tempIdx),2);
        stdData{m}=std(rnaSeqMat(:,tempIdx),[],2);
        meanNames{m}=[strains{i} ' ' conditions{j}];
        m=m+1;
        
    end
    
end

%take mean of untreated parents
vTpm=mean([meanData{1} meanData{2}],2);
vTpm(isnan(vTpm))=1;    %set missing to 1 TPM

[sortedTpm,sortIdx]=sort(vTpm,'ascend');
sortedTpmGenes=geneNames(sortIdx);


%expressionFlag='_topQuintile';
%expressionFlag='_bottomQuintile';
expressionFlag='_allGenes';

doSim=1;


tempLength=floor(length(sortedTpmGenes)/5);
%expressionGenesToUse=sortedTpmGenes((end-tempLength):end);
expressionGenesToUse=sortedTpmGenes(1:tempLength);

idxToUse=ismember(gene,expressionGenesToUse);
%idxToUse=logical(ones(length(gene),1));






% queryGene=[];
% secondary=cell(length(gene),1);
% sasa=nan(length(gene),1);
% normSasa=nan(length(gene),1);
% refCodon=nan(length(gene),1);
% altCodon=nan(length(gene),1);
% phi=nan(length(gene),1);
% psi=nan(length(gene),1);
% posInOrf=nan(length(gene),1);
% posInRun=nan(length(gene),1);
% lengthOfRun=nan(length(gene),1);
% neighbors=nan(length(gene),1);
% pDist=nan(length(gene),1);
% for i=1:length(gene)
%     
%     queryGene=gene{i};
%       
%     toGet1=['mutationOutput/' queryGene 'dsspTable.mat'];
%     toGet2=['mutationOutput/' queryGene 'neighborTable.mat'];
%     
%     if sum(ismember(systematicNamePdb,queryGene))>0
%         pdbQuery=pdbId(ismember(systematicNamePdb,queryGene));
%         pdbQuery=pdbQuery{1};   %in case of multiple mapping
%         toGet3=['mutationOutput/' pdbQuery 'pDist.mat'];
%     else
%         toGet3='NA';
%     end
% 
%     if exist(toGet1)&&exist(toGet2)
% 
%         load(toGet1)
%         dsspTable=outputTable;
%         
%         load(toGet2)
%         neighborTable=outputTable;
%         
%         if exist(toGet3)
%             load(toGet3)
%             pDistTable=outputTable;
%         end
% 
%         %grab secondary structure and sasa
%         residueQuery=str2double(proteinEncoded{i}(6:(end-3)));
%         baseQuery=str2double(dnaEncoded{i}(3:(end-3)));
%         
%         %get ref and alt codons
%         baseQueryRange=(3*residueQuery-2):(3*residueQuery);
%         baseToChange=find(ismember(baseQueryRange,baseQuery));
%         
%         if max(baseQueryRange)<=length(dnaSequence)
%             refBase=dnaEncoded{i}(end-2);
%             altBase=dnaEncoded{i}(end);
% 
%             tempRefCodon=dnaSequence(baseQueryRange);
%             refCodon(i)=codonTypes({tempRefCodon});
% 
%             tempAltCodon=tempRefCodon;
%             tempAltCodon(baseToChange)=altBase;
%             altCodon(i)=codonTypes({tempAltCodon});
%         end
% 
%         if residueQuery<=height(dsspTable)
%             secondary{i}=dsspTable.secondary{residueQuery};
%             sasa(i)=dsspTable.sasa(residueQuery);
%             normSasa(i)=sasa(i)/mean(dsspTable.sasa,'omitnan');
%             phi(i)=dsspTable.phi(residueQuery);
%             psi(i)=dsspTable.psi(residueQuery);
%             posInOrf(i)=dsspTable.posInOrf(residueQuery);
%             posInRun(i)=dsspTable.posInRun(residueQuery);
%             lengthOfRun(i)=dsspTable.lengthOfRun(residueQuery);
%             neighbors(i)=neighborTable.neighbors(residueQuery);
%         end
%         
%         if exist(toGet3)&&(residueQuery<length(pDistTable.pDist))
%             pDist(i)=pDistTable.pDist(residueQuery);
%         else
%             pDist(i)=nan;
%         end
%         
%     end
%     
%     if mod(i,1000)==0
%         i
%     end
%     
% end
% 
% %save output to avoid rerunning
% save('secondary1k.mat','secondary')
% save('sasa1k.mat','sasa')
% save('normSasa1k.mat','normSasa')
% save('phi1k.mat','phi')
% save('psi1k.mat','psi')
% save('posInOrf1k.mat','posInOrf')
% save('posInRun1k.mat','posInRun')
% save('lengthOfRun1k.mat','lengthOfRun')
% save('refCodon1k.mat','refCodon')
% save('altCodon1k.mat','altCodon')
% save('neighbors1k.mat','neighbors')
% save('pDist1k.mat','pDist')
% 






load('secondary1k.mat','secondary')
load('sasa1k.mat','sasa')
load('normSasa1k.mat','normSasa')
load('phi1k.mat','phi')
load('psi1k.mat','psi')
load('posInOrf1k.mat','posInOrf')
load('posInRun1k.mat','posInRun')
load('lengthOfRun1k.mat','lengthOfRun')
load('refCodon1k.mat','refCodon')
load('altCodon1k.mat','altCodon')
load('neighbors1k.mat','neighbors')
load('pDist1k.mat','pDist')

load('simulationStructureData.mat')

%clear out some unnecessary variables
clear chr dnaEncoded dnaSequence dsspTable inputGenome %gene ...
clear mutationTable outputTable pdbDictionary pdbId

%organize data to plot systematically
propertyLabels={'accessible surface area','normalized ASA','phi','psi','position in ORF',...
    'position in structure run','length of structure run',...
    'relative position in structure run','ref nTE','delta nTE','|delta nTE|',...
    'neighbors','dist to Psite'};

%xProperties=[1:2 5:8];
yProperties=1:length(propertyLabels);
xProperties=yProperties;

yLim1=[0 0 -180 -180 0 0 0 0 0 -0.5 0 0 0];
%yLim2=[150 2 180 180 1 50 50 1 0.5 0.5 0.5 25 50];
yLim2=[250 2 180 180 1 50 50 1 0.5 0.5 0.5 40 100];
%make array of properties for 1k genomes and simulations
%also secondary structure for sorting and isTS
%for sim, do mis and syn separately as well as combined
%do combined and then use vMis and vSyn to separate?

%1k first
type=type(idxToUse);
vMis1k=ismember(type,'missense_variant');
vSyn1k=ismember(type,'synonymous_variant');

clear type

secondary=secondary(idxToUse);
vHasStruct=~cellfun(@isempty,secondary);

[~,structureMis1k]=structureTypes(secondary(logical(vMis1k.*vHasStruct)));  %for categorizing

[~,structureSyn1k]=structureTypes(secondary(logical(vSyn1k.*vHasStruct)));  %for categorizing

clear secondary

sasa=sasa(idxToUse);
normSasa=normSasa(idxToUse);
phi=phi(idxToUse);
psi=psi(idxToUse);
posInOrf=posInOrf(idxToUse);
posInRun=posInRun(idxToUse);
lengthOfRun=lengthOfRun(idxToUse);
refCodon=refCodon(idxToUse);
altCodon=altCodon(idxToUse);
neighbors=neighbors(idxToUse);
pDist=pDist(idxToUse);

properties1k{1,1}=sasa(logical(vMis1k.*vHasStruct));
properties1k{1,2}=normSasa(logical(vMis1k.*vHasStruct));
properties1k{1,3}=phi(logical(vMis1k.*vHasStruct));
properties1k{1,4}=psi(logical(vMis1k.*vHasStruct));
properties1k{1,5}=posInOrf(logical(vMis1k.*vHasStruct));
properties1k{1,6}=posInRun(logical(vMis1k.*vHasStruct));
properties1k{1,7}=lengthOfRun(logical(vMis1k.*vHasStruct));
properties1k{1,8}=posInRun(logical(vMis1k.*vHasStruct))./lengthOfRun(logical(vMis1k.*vHasStruct));
[~,tempRefnTE]=codonFrequency(refCodon(logical(vMis1k.*vHasStruct)));
[~,tempAltnTE]=codonFrequency(altCodon(logical(vMis1k.*vHasStruct)));
properties1k{1,9}=tempRefnTE;
properties1k{1,10}=tempAltnTE-tempRefnTE;
properties1k{1,11}=abs(tempAltnTE-tempRefnTE);
properties1k{1,12}=neighbors(logical(vMis1k.*vHasStruct));
properties1k{1,13}=pDist(logical(vMis1k.*vHasStruct));


properties1k{2,1}=sasa(logical(vSyn1k.*vHasStruct));
properties1k{2,2}=normSasa(logical(vSyn1k.*vHasStruct));
properties1k{2,3}=phi(logical(vSyn1k.*vHasStruct));
properties1k{2,4}=psi(logical(vSyn1k.*vHasStruct));
properties1k{2,5}=posInOrf(logical(vSyn1k.*vHasStruct));
properties1k{2,6}=posInRun(logical(vSyn1k.*vHasStruct));
properties1k{2,7}=lengthOfRun(logical(vSyn1k.*vHasStruct));
properties1k{2,8}=posInRun(logical(vSyn1k.*vHasStruct))./lengthOfRun(logical(vSyn1k.*vHasStruct));
[~,tempRefnTE]=codonFrequency(refCodon(logical(vSyn1k.*vHasStruct)));
[~,tempAltnTE]=codonFrequency(altCodon(logical(vSyn1k.*vHasStruct)));
properties1k{2,9}=tempRefnTE;
properties1k{2,10}=tempAltnTE-tempRefnTE;
properties1k{2,11}=abs(tempAltnTE-tempRefnTE);
properties1k{2,12}=neighbors(logical(vSyn1k.*vHasStruct));
properties1k{2,13}=pDist(logical(vSyn1k.*vHasStruct));


%allele frequencies
af=af(idxToUse);
misAf=af(logical(vMis1k.*vHasStruct));
synAf=af(logical(vSyn1k.*vHasStruct));


clear af sasa phi psi posInOrf posInRun lengthOfRun neighbors pDist %altCodon refCodon

%keep ref and alt from 1K to do per-codon analysis




structureMisSim=misSecondary;   %for categorizing

propertiesSim{1,1}=misSasa;
propertiesSim{1,2}=misNormSasa;
propertiesSim{1,3}=misPhi;
propertiesSim{1,4}=misPsi;
propertiesSim{1,5}=misPosInOrf;
propertiesSim{1,6}=misPosInRun;
propertiesSim{1,7}=misLengthOfRun;
propertiesSim{1,8}=misPosInRun./misLengthOfRun;
[~,tempRefnTE]=codonFrequency(misInitialCodon);
[~,tempAltnTE]=codonFrequency(misNewCodon);
propertiesSim{1,9}=tempRefnTE;
propertiesSim{1,10}=tempAltnTE-tempRefnTE;
propertiesSim{1,11}=abs(tempAltnTE-tempRefnTE);
propertiesSim{1,12}=misNeighbors;
propertiesSim{1,13}=misPdist;

structureSynSim=synSecondary;   %for categorizing

propertiesSim{2,1}=synSasa;
propertiesSim{2,2}=synNormSasa;
propertiesSim{2,3}=synPhi;
propertiesSim{2,4}=synPsi;
propertiesSim{2,5}=synPosInOrf;
propertiesSim{2,6}=synPosInRun;
propertiesSim{2,7}=synLengthOfRun;
propertiesSim{2,8}=synPosInRun./synLengthOfRun;
[~,tempRefnTE]=codonFrequency(synInitialCodon);
[~,tempAltnTE]=codonFrequency(synNewCodon);
propertiesSim{2,9}=tempRefnTE;
propertiesSim{2,10}=tempAltnTE-tempRefnTE;
propertiesSim{2,11}=abs(tempAltnTE-tempRefnTE);
propertiesSim{2,12}=synNeighbors;
propertiesSim{2,13}=synPdist;

clear misSasa misPhi misPsi misPosInOrf misPosInRun misLengthOfRun misNeighbors misInitialCodon misNewCodon 
clear synSasa synPhi synPsi synPosInOrf synPosInRun synLengthOfRun synNeighbors %synInitialCodon synNewCodon 

toc

afThresh=1e-3;




%codon array from codonTypes
bases={'A','C','G','T'};

codons=cell(64,1);
m=1;
for i=1:length(bases)
    
    for j=1:length(bases)
        
        for k=1:length(bases)
            
            codons{m}=[bases{i} bases{j} bases{k}];
            m=m+1;
            
        end
        
    end
    
end

structToUse=[1:2 5:8];


%simplified figures for paper
%leave out Ts/Tv comparison -- just for supp
%use median normalization and error bars?
for i=[1 12 13]%1:length(propertyLabels)

    tic
    
    if doSim
    %gross distribution (before considering stucture)
    figure('units','normalized','outerposition',[0 0 1 1])

    
    %mis then syn; merge to make all
    v1=properties1k{1,i};
    v2=properties1k{2,i};
    v3=[v1;v2];
    
    
    v4=propertiesSim{1,i};
    v5=propertiesSim{2,i};
    v6=[v4;v5];
    
    clear toPlot
    %all
    toPlot{1}=v3;
    toPlot{2}=v6;
    %toPlot{3}=v6(logical([misTs;synTs]));
    %toPlot{4}=v6(~logical([misTs;synTs]));
    
    %tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    tempLabels={'observed 1k','simulated'};
    subplot(2,6,1)
    easyBox(toPlot)
    
    
    
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    %missense
    toPlot{1}=v1;
    toPlot{2}=v4;
    %toPlot{3}=v4(logical(misTs));
    %toPlot{4}=v4(~logical(misTs));
    
    %tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    tempLabels={'observed 1k','simulated'};
    subplot(2,6,2)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    %histograms
    subplot(2,4,3)
    hold on
    histogram(toPlot{1},0:(ceil(yLim2(i)/20)):(yLim2(i)),'Normalization','Probability')
    histogram(toPlot{2},0:(ceil(yLim2(i)/20)):(yLim2(i)),'Normalization','Probability')
    xlim([0 yLim2(i)])
    %xlim([0 500])
    axis square
    title('missense')
    
    
    
    %synonymous
    toPlot{1}=v2;
    toPlot{2}=v5;
    %toPlot{3}=v5(logical(synTs));
    %toPlot{4}=v5(~logical(synTs));
    
    %tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    tempLabels={'observed 1k','simulated'};
    subplot(2,6,3)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    subplot(2,4,4)
    hold on
    histogram(toPlot{1},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    histogram(toPlot{2},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    xlim([0 yLim2(i)])
    axis square
    title('synonymous')
    
    
    
%     %break down by secondary structure for all;mis;syn
%     clear toPlot
%     clear tempLabels
%     m=1;
%     for j=1:length(structureLabels(structToUse))
%         
%         %observed
%         tempIdx1=[structureMis1k';structureSyn1k']==j;
%         
%         %sim
%         tempIdx2=[misSecondary;synSecondary]==j;
%         
%         toPlot(m)=median(v3(tempIdx1),'omitnan')/median(v6(tempIdx2),'omitnan');
%         tempLabels{m}=['obs./sim. ' structureLabels{structToUse(j)}];
%         m=m+1;
%         
%     end
%     
%     subplot(2,2,2)
%     hold on
%     %easyBox(toPlot)
%     bar(toPlot,'BaseValue',1)
%     xticks(1:length(tempLabels))
%     xtickangle(45)
%     xticklabels(tempLabels)
%     ylabel(propertyLabels{i})
%     title('missense and synonymous')
%     %ylim([yLim1(i) yLim2(i)])
%     xlim([0.5 length(tempLabels)+0.5])
%     ylim([1/2 2])
%     set(gca,'YScale','log')
    
    
    
    
    clear toPlot tempLabels pVal
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        
        %sim
        tempIdx2=misSecondary==j;
        
        toPlot(m)=median(v1(tempIdx1),'omitnan')/median(v4(tempIdx2),'omitnan');
        tempLabels{m}=['obs./sim. ' structureLabels{structToUse(j)}];
        [pVal(m) h]=ranksum(v1(tempIdx1),v4(tempIdx2));
        m=m+1;
        
        
    end
    
    subplot(2,2,3)
    hold on
    %easyBox(toPlot)
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    %ylim([yLim1(i) yLim2(i)])
    xlim([0.5 length(tempLabels)+0.5])
    ylim([1/2 2])
    set(gca,'YScale','log')
    for j=1:length(pVal)
        text(j,1.9,num2str(pVal(j)))
    end
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureSyn1k'==j;
        
        %sim
        tempIdx2=synSecondary==j;
        
        toPlot(m)=median(v2(tempIdx1),'omitnan')/median(v5(tempIdx2),'omitnan');
        tempLabels{m}=['obs./sim. ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,2,4)
    hold on
    %easyBox(toPlot)
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    %ylim([yLim1(i) yLim2(i)])
    xlim([0.5 length(tempLabels)+0.5])
    ylim([1/2 2])
    set(gca,'YScale','log')
    
    
    set(gcf,'PaperPositionMode','auto')
    print(['figures/structureSelection_' propertyLabels{i} '_fig_' num2str(figureCounter)],'-dsvg','-r0')
    print(['figures/structureSelection_' propertyLabels{i} '_fig_' num2str(figureCounter)],'-djpeg','-r0')
    figureCounter=figureCounter+1;
    
    
    
    %close
    
    end
    
    
    %rare vs common odds ratio plots?
    clear toPlot
    figure('units','normalized','outerposition',[0 0 1 1])

    %mis then syn; merge to make all
    v1=properties1k{1,i};
    v2=properties1k{2,i};
    v3=[v1;v2];
    
    misIdx=misAf<afThresh;
    synIdx=synAf<afThresh;
    
    toPlot{1}=v3(logical([misIdx;synIdx]));
    toPlot{2}=v3(~logical([misIdx;synIdx]));
    
    tempLabels={'rare in 1k','common in 1k'};
    subplot(2,6,1)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    toPlot{1}=v1(misIdx);
    toPlot{2}=v1(~misIdx);
    
    subplot(2,6,2)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    %also histogram
    subplot(2,4,3)
    hold on
    histogram(toPlot{1},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    histogram(toPlot{2},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    xlim([0 yLim2(i)])
    axis square
    title('missense')
    
    
    
    toPlot{1}=v2(synIdx);
    toPlot{2}=v2(~synIdx);
    
    subplot(2,6,3)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    %subset by structure
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=[structureMis1k';structureSyn1k']==j;
        tempIdx2=logical([misIdx;synIdx]);
        
        toPlot(m)=median(v3(logical(tempIdx1.*tempIdx2)),'omitnan')/...
            median(v3(logical(tempIdx1.*~tempIdx2)),'omitnan');
        tempLabels{m}=['rare/common' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,4)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    xlim([0.5 length(tempLabels)+0.5])
    ylim([1/2 2])
    set(gca,'YScale','log')
    title('missense and synonymous')
    
    
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        tempIdx2=misIdx;
        
        toPlot(m)=median(v1(logical(tempIdx1.*tempIdx2)),'omitnan')/...
            median(v1(logical(tempIdx1.*~tempIdx2)),'omitnan');
        tempLabels{m}=['rare/common ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,5)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    xlim([0.5 length(tempLabels)+0.5])
    ylim([2/3 3/2])
    set(gca,'YScale','log')
    title('missense')
    %pVals
    
    
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureSyn1k'==j;
        tempIdx2=synIdx;
        
        toPlot(m)=median(v2(logical(tempIdx1.*tempIdx2)),'omitnan')/...
            median(v2(logical(tempIdx1.*~tempIdx2)),'omitnan');
        tempLabels{m}=['rare/common ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,6)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    xlim([0.5 length(tempLabels)+0.5])
    ylim([2/3 3/2])
    set(gca,'YScale','log')
    title('synonymous')
    
    
        
    
    
    set(gcf,'PaperPositionMode','auto')
    print(['figures/structureSelection_' propertyLabels{i} expressionFlag '_fig_' num2str(figureCounter)],'-dsvg','-r0')
    print(['figures/structureSelection_' propertyLabels{i} expressionFlag '_fig_' num2str(figureCounter)],'-djpeg','-r0')
    figureCounter=figureCounter+1;
    
    
    %close
    
    
    toc
    
    
    
end


jyntbsre

for i=1:length(propertyLabels)

    tic
    
    if doSim
    %gross distribution (before considering stucture)
    figure('units','normalized','outerposition',[0 0 1 1])

    
    %mis then syn; merge to make all
    v1=properties1k{1,i};
    v2=properties1k{2,i};
    v3=[v1;v2];
    
    
    v4=propertiesSim{1,i};
    v5=propertiesSim{2,i};
    v6=[v4;v5];
    
    clear toPlot
    %all
    toPlot{1}=v3;
    toPlot{2}=v6;
    toPlot{3}=v6(logical([misTs;synTs]));
    toPlot{4}=v6(~logical([misTs;synTs]));
    
    tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    subplot(2,6,1)
    easyBox(toPlot)
    
    
    
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    %missense
    toPlot{1}=v1;
    toPlot{2}=v4;
    toPlot{3}=v4(logical(misTs));
    toPlot{4}=v4(~logical(misTs));
    
    tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    subplot(2,6,2)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    %synonymous
    toPlot{1}=v2;
    toPlot{2}=v5;
    toPlot{3}=v5(logical(synTs));
    toPlot{4}=v5(~logical(synTs));
    
    tempLabels={'observed 1k','simulated','sim. only Ts','sim. only Tv'};
    subplot(2,6,3)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    %break down by secondary structure for all;mis;syn
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=[structureMis1k';structureSyn1k']==j;
        toPlot{m}=v3(tempIdx1);
        tempLabels{m}=['observed 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim
        tempIdx2=[misSecondary;synSecondary]==j;
        toPlot{m}=v6(tempIdx2);
        tempLabels{m}=['simulated ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Ts
        tempIdx3=[misTs;synTs]==1;
        toPlot{m}=v6(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Ts ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Tv
        tempIdx3=[misTs;synTs]~=1;
        toPlot{m}=v6(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Tv ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,2,2)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        toPlot{m}=v1(tempIdx1);
        tempLabels{m}=['observed 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim
        tempIdx2=misSecondary==j;
        toPlot{m}=v4(tempIdx2);
        tempLabels{m}=['simulated ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Ts
        tempIdx3=misTs==1;
        toPlot{m}=v4(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Ts ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Tv
        tempIdx3=misTs~=1;
        toPlot{m}=v4(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Tv ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,2,3)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureSyn1k'==j;
        toPlot{m}=v2(tempIdx1);
        tempLabels{m}=['observed 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim
        tempIdx2=synSecondary==j;
        toPlot{m}=v5(tempIdx2);
        tempLabels{m}=['simulated ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Ts
        tempIdx3=synTs==1;
        toPlot{m}=v5(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Ts ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %sim Tv
        tempIdx3=synTs~=1;
        toPlot{m}=v5(logical(tempIdx2.*tempIdx3));
        tempLabels{m}=['simulated only Tv ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,2,4)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    
    set(gcf,'PaperPositionMode','auto')
    print(['figures/structureSelection_' propertyLabels{i} '_fig_' num2str(figureCounter)],'-dsvg','-r0')
    print(['figures/structureSelection_' propertyLabels{i} '_fig_' num2str(figureCounter)],'-djpeg','-r0')
    figureCounter=figureCounter+1;
    
    close
    
    end
    
    
    %also plot median ratios
    
    ntbrawvra
    

    %rare vs common odds ratio plots?
    clear toPlot
    figure('units','normalized','outerposition',[0 0 1 1])

    %mis then syn; merge to make all
    v1=properties1k{1,i};
    v2=properties1k{2,i};
    v3=[v1;v2];
    
    misIdx=misAf<afThresh;
    synIdx=synAf<afThresh;
    
    toPlot{1}=v3(logical([misIdx;synIdx]));
    toPlot{2}=v3(~logical([misIdx;synIdx]));
    
    tempLabels={'rare in 1k','common in 1k'};
    subplot(2,8,1)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    
    toPlot{1}=v1(misIdx);
    toPlot{2}=v1(~misIdx);
    
    subplot(2,8,2)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    toPlot{1}=v2(synIdx);
    toPlot{2}=v2(~synIdx);
    
    subplot(2,8,3)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    n=1;
    for ii=1:length(toPlot)
        for jj=(ii+1):length(toPlot)
            if (sum(~isnan(toPlot{ii}))>0)&&(sum(~isnan(toPlot{jj}))>0)
                [p h]=ranksum(toPlot{ii},toPlot{jj});
                plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
                text((ii+jj)/2,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
                n=n+1;
            end
        end
    end
    
    %subset by structure
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=[structureMis1k';structureSyn1k']==j;
        tempIdx2=logical([misIdx;synIdx]);
        toPlot{m}=v3(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['rare in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %observed
        tempIdx2=~logical([misIdx;synIdx]);
        toPlot{m}=v3(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['common in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,4)
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense and synonymous')
    ylim([yLim1(i) yLim2(i)])
    hold on
    
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        tempIdx2=misIdx;
        toPlot{m}=v1(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['rare in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %observed
        tempIdx2=~misIdx;
        toPlot{m}=v1(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['common in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,5)
    hold on
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('missense')
    ylim([yLim1(i) yLim2(i)])
    %pVals
    n=1;
    for ii=1:2:length(toPlot)
        [p h]=ranksum(toPlot{ii},toPlot{ii+1});
        plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
        text(ii+0.5,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
        n=n+1;
    end
    
    
    
    clear toPlot
    clear tempLabels
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureSyn1k'==j;
        tempIdx2=synIdx;
        toPlot{m}=v2(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['rare in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
        %observed
        tempIdx2=~synIdx;
        toPlot{m}=v2(logical(tempIdx1.*tempIdx2));
        tempLabels{m}=['common in 1k ' structureLabels{structToUse(j)}];
        m=m+1;
        
    end
    
    subplot(2,3,6)
    hold on
    easyBox(toPlot)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel(propertyLabels{i})
    title('synonymous')
    ylim([yLim1(i) yLim2(i)])
    %pVals
    n=1;
    for ii=1:2:length(toPlot)
        [p h]=ranksum(toPlot{ii},toPlot{ii+1});
        plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
        text(ii+0.5,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
        n=n+1;
    end
    
    
    
    
    
    set(gcf,'PaperPositionMode','auto')
    print(['figures/structureSelection_' propertyLabels{i} expressionFlag '_fig_' num2str(figureCounter)],'-dsvg','-r0')
    print(['figures/structureSelection_' propertyLabels{i} expressionFlag '_fig_' num2str(figureCounter)],'-djpeg','-r0')
    figureCounter=figureCounter+1;
    
    
    close
    
    
    toc
    
    
    rrbavew
    
    
end


%codon combos; just do select properties

%also do this analysis for sim [codon pairs should reduce Ts/Tv issues]

%not enough combos in 1K;skip

% for l=1:64  %ref
%     
%     for k=1:64  %alt
%         
%         
%             refAltIdx=logical((refCodon(logical(vSyn1k.*vHasStruct))==l).*...
%                 (altCodon(logical(vSyn1k.*vHasStruct))==k));
%             
%             simRefAltIdx=logical((synInitialCodon==l).*...
%                 (synNewCodon==k));
%             
%             
%             
%             if sum(refAltIdx)>100
%             
%                 synIdx=synAf<afThresh;
%             
%                 figure('units','normalized','outerposition',[0 0 1 1])
% 
%                 o=1;
% 
%                 for i=[1 12 5]
%                     
%                     %mis then syn; merge to make all
%                     v1=properties1k{1,i};
%                     v2=properties1k{2,i};
%                     v3=[v1;v2];
%     
%                     v4=propertiesSim{1,i};
%                     v5=propertiesSim{2,i};
%                     v6=[v4;v5];
% 
%                     
%                     clear toPlot
%                     clear tempLabels
%                     m=1;
%                     for j=1:length(structureLabels(structToUse))
% 
%                         %observed
%                         tempIdx1=structureSyn1k'==j;
%                         tempIdx2=logical(synIdx.*refAltIdx);
%                         toPlot{m}=v2(logical(tempIdx1.*tempIdx2));
%                         tempLabels{m}=['observed 1k ' structureLabels{j}];
%                         m=m+1;
% 
%                         %sim
%                         tempIdx1=synSecondary==j;
%                         tempIdx2=simRefAltIdx;
%                         toPlot{m}=v5(logical(tempIdx1.*tempIdx2));
%                         tempLabels{m}=['simulated ' structureLabels{j}];
%                         m=m+1;
% 
%                         %sim Ts
%                         tempIdx3=synTs==1;
%                         toPlot{m}=v5(logical(tempIdx1.*tempIdx2.*tempIdx3));
%                         tempLabels{m}=['simulated only Ts ' structureLabels{structToUse(j)}];
%                         m=m+1;
% 
%                         %sim Tv
%                         tempIdx3=synTs~=1;
%                         toPlot{m}=v5(logical(tempIdx1.*tempIdx2.*tempIdx3));
%                         tempLabels{m}=['simulated only Tv ' structureLabels{structToUse(j)}];
%                         m=m+1;
% 
%                     end
% 
%                     subplot(2,3,o)
%                     easyBox(toPlot)
%                     xticks(1:length(tempLabels))
%                     xtickangle(45)
%                     xticklabels(tempLabels)
%                     ylabel(propertyLabels{i})
%                     title('synonymous')
%                     ylim([yLim1(i) yLim2(i)])
%     
% 
%                     
% 
% 
%                     clear toPlot
%                     clear tempLabels
%                     m=1;
%                     for j=1:length(structureLabels(structToUse))
% 
%                         %observed
%                         tempIdx1=structureSyn1k'==j;
%                         tempIdx2=logical(synIdx.*refAltIdx);
%                         toPlot{m}=v2(logical(tempIdx1.*tempIdx2));
%                         tempLabels{m}=['rare in 1k ' structureLabels{structToUse(j)}];
%                         m=m+1;
% 
%                         %observed
%                         tempIdx2=logical(~synIdx.*refAltIdx);
%                         toPlot{m}=v2(logical(tempIdx1.*tempIdx2));
%                         tempLabels{m}=['common in 1k ' structureLabels{structToUse(j)}];
%                         m=m+1;
% 
%                     end
% 
%                     subplot(2,3,o+3)
%                     o=o+1;
%                     hold on
%                     easyBox(toPlot)
%                     xticks(1:length(tempLabels))
%                     xtickangle(45)
%                     xticklabels(tempLabels)
%                     ylabel(propertyLabels{i})
%                     title([codons{l} ' ' codons{k} ' synonymous'])
%                     ylim([yLim1(i) yLim2(i)])
%                     %pVals
%                     n=1;
%                     for ii=1:2:length(toPlot)
%                         jj=ii+1;
%                         if (length(toPlot{ii})>0)&&(length(toPlot{ii+1})>0)
%                             [p h]=ranksum(toPlot{ii},toPlot{ii+1});
%                             plot([ii jj],[(0.6+(0.05*(n-1)))*yLim2(i) (0.6+(0.05*(n-1)))*yLim2(i)],'-k')
%                             text(ii+0.5,(0.6+(0.05*(n-1))+0.025)*yLim2(i),num2str(p))
%                         end
%                         n=n+1;
%                     end
% 
% 
%                 end
% 
%                 
%                 set(gcf,'PaperPositionMode','auto')
%                 print(['figures/structureSelection_' codons{l} '_' codons{k} '_fig_' num2str(figureCounter)],'-dsvg','-r0')
%                 print(['figures/structureSelection_' codons{l} '_' codons{k} '_fig_' num2str(figureCounter)],'-djpeg','-r0')
%                 figureCounter=figureCounter+1;
% 
%                 close
%                 
%             end
%         
%     end
%     
% end





%figureCounter=23;




%make some 2-D plots by position
nBins=10;
for i=1:length(xProperties)
    
    %mis then syn; merge to make all
    x1=properties1k{1,xProperties(i)};
    x2=properties1k{2,xProperties(i)};
    x3=[x1;x2];
    
    x4=propertiesSim{1,xProperties(i)};
    x5=propertiesSim{2,xProperties(i)};
    x6=[x4;x5];

    
    for j=1:length(yProperties)

        tic

        %gross distribution (before considering stucture)
        figure('units','normalized','outerposition',[0 0 1 1])


        y1=properties1k{1,yProperties(j)};
        y2=properties1k{2,yProperties(j)};
        y3=[y1;y2];
        
        y4=propertiesSim{1,yProperties(j)};
        y5=propertiesSim{2,yProperties(j)};
        y6=[y4;y5];

        %add only Ts and only Tv
        
        [tempMedian1,~]=getDeciles(x3,y3,nBins);
        [tempMedian2,~]=getDeciles(x6,y6,nBins);
        [tempMedian3,~]=getDeciles(x6(logical([misTs;synTs])),...
            y6(logical([misTs;synTs])),nBins);
        [tempMedian4,~]=getDeciles(x6(~logical([misTs;synTs])),...
            y6(~logical([misTs;synTs])),nBins);
        tempMedian5=(3*tempMedian3+tempMedian4)/4;
        
        %also subset by structure for mis and syn (just do top plot?)
        
        subplot(2,3,1)
        hold on
        plot(1:length(tempMedian1),tempMedian1,'Color',blue)
        plot(1:length(tempMedian2),tempMedian2,'Color',orange)
        plot(1:length(tempMedian3),tempMedian3,'Color',grey)
        plot(1:length(tempMedian4),tempMedian4,'Color','black')
        plot(1:length(tempMedian5),tempMedian5,':','Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([yLim1(yProperties(j)) yLim2(yProperties(j))])
        title('missense and synonymous')
        xlim([1 nBins])
        axis square
        legend({'observed','simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        %relative
        subplot(2,3,4)
        hold on
        plot(1:length(tempMedian1),tempMedian1./tempMedian2,'Color',blue)
        plot(1:length(tempMedian1),tempMedian1./tempMedian3,'Color',grey)
        plot(1:length(tempMedian1),tempMedian1./tempMedian4,':','Color','black')
        plot(1:length(tempMedian1),tempMedian1./tempMedian5,'Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([0.5 2])
        set(gca,'YScale','log')
        title('missense and synonymous')
        xlim([1 nBins])
        axis square
        plot(xlim,[1 1],':k')
        legend({'simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        
        
        
        [tempMedian1,~]=getDeciles(x1,y1,nBins);
        [tempMedian2,~]=getDeciles(x4,y4,nBins);
        [tempMedian3,~]=getDeciles(x4(logical(misTs)),y4(logical(misTs)),nBins);
        [tempMedian4,~]=getDeciles(x4(~logical(misTs)),y4(~logical(misTs)),nBins);
        tempMedian5=(3*tempMedian3+tempMedian4)/4;
        
        subplot(2,3,2)
        hold on
        plot(1:length(tempMedian1),tempMedian1,'Color',blue)
        plot(1:length(tempMedian2),tempMedian2,'Color',orange)
        plot(1:length(tempMedian3),tempMedian3,'Color',grey)
        plot(1:length(tempMedian4),tempMedian4,':','Color','black')
        plot(1:length(tempMedian5),tempMedian5,'Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([yLim1(yProperties(j)) yLim2(yProperties(j))])
        title('missense')
        xlim([1 nBins])
        axis square
        legend({'observed','simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        %relative
        subplot(2,3,5)
        hold on
        plot(1:length(tempMedian1),tempMedian1./tempMedian2,'Color',blue)
        plot(1:length(tempMedian1),tempMedian1./tempMedian3,'Color',grey)
        plot(1:length(tempMedian1),tempMedian1./tempMedian4,':','Color','black')
        plot(1:length(tempMedian1),tempMedian1./tempMedian5,'Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([0.5 2])
        set(gca,'YScale','log')
        title('missense')
        xlim([1 nBins])
        axis square
        plot(xlim,[1 1],':k')
        legend({'simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        
        
        
        [tempMedian1,~]=getDeciles(x2,y2,nBins);
        [tempMedian2,~]=getDeciles(x5,y5,nBins);
        [tempMedian3,~]=getDeciles(x5(logical(synTs)),y5(logical(synTs)),nBins);
        [tempMedian4,~]=getDeciles(x5(~logical(synTs)),y5(~logical(synTs)),nBins);
        tempMedian5=(3*tempMedian3+tempMedian4)/4;
        
        subplot(2,3,3)
        hold on
        plot(1:length(tempMedian1),tempMedian1,'Color',blue)
        plot(1:length(tempMedian2),tempMedian2,'Color',orange)
        plot(1:length(tempMedian3),tempMedian3,'Color',grey)
        plot(1:length(tempMedian4),tempMedian4,':','Color','black')
        plot(1:length(tempMedian5),tempMedian5,'Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([yLim1(yProperties(j)) yLim2(yProperties(j))])
        title('synonymous')
        xlim([1 nBins])
        axis square
        legend({'observed','simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        %relative
        subplot(2,3,6)
        hold on
        plot(1:length(tempMedian1),tempMedian1./tempMedian2,'Color',blue)
        plot(1:length(tempMedian1),tempMedian1./tempMedian3,'Color',grey)
        plot(1:length(tempMedian1),tempMedian1./tempMedian4,':','Color','black')
        plot(1:length(tempMedian1),tempMedian1./tempMedian5,'Color','black')
        xlabel(propertyLabels{xProperties(i)})
        ylabel(propertyLabels{yProperties(j)})
        ylim([0.5 2])
        set(gca,'YScale','log')
        title('synonymous')
        xlim([1 nBins])
        axis square
        plot(xlim,[1 1],':k')
        legend({'simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
        
        
        set(gcf,'PaperPositionMode','auto')
        print(['figures/structureSelection_' propertyLabels{xProperties(i)} '_' propertyLabels{yProperties(j)} '_fig_' num2str(figureCounter)],'-dsvg','-r0')
        print(['figures/structureSelection_' propertyLabels{xProperties(i)} '_' propertyLabels{yProperties(j)} '_fig_' num2str(figureCounter)],'-djpeg','-r0')
        figureCounter=figureCounter+1;
        
        close

        figure('units','normalized','outerposition',[0 0 1 1])
        %missense by structure
        for l=1:length(structureLabels)
            
            tempIdx1=structureMis1k==l;
            tempIdx2=structureMisSim==l;
            tempIdx3=logical(misTs);
            
            [tempMedian1,~]=getDeciles(x1(tempIdx1),y1(tempIdx1),nBins);
            [tempMedian2,~]=getDeciles(x4(tempIdx2),y4(tempIdx2),nBins);
            [tempMedian3,~]=getDeciles(x4(logical(tempIdx2.*tempIdx3)),y4(logical(tempIdx2.*tempIdx3)),nBins);
            [tempMedian4,~]=getDeciles(x4(logical(tempIdx2.*~tempIdx3)),y4(logical(tempIdx2.*~tempIdx3)),nBins);
            tempMedian5=(3*tempMedian3+tempMedian4)/4;
            
            subplot(4,4,l)
            hold on
            plot(1:length(tempMedian1),tempMedian1,'Color',blue)
            plot(1:length(tempMedian2),tempMedian2,'Color',orange)
            plot(1:length(tempMedian3),tempMedian3,'Color',grey)
            plot(1:length(tempMedian4),tempMedian4,':','Color','black')
            plot(1:length(tempMedian5),tempMedian5,'Color','black')
            xlabel(propertyLabels{xProperties(i)})
            ylabel(propertyLabels{yProperties(j)})
            ylim([yLim1(yProperties(j)) yLim2(yProperties(j))])
            title(['missense ' structureLabels{l}])
            xlim([1 nBins])
            axis square
            %legend({'observed','simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
            
        end

        %synonymous by structure
        for l=1:length(structureLabels)
            
            tempIdx1=structureSyn1k==l;
            tempIdx2=structureSynSim==l;
            tempIdx3=logical(synTs);
            
            [tempMedian1,~]=getDeciles(x2(tempIdx1),y2(tempIdx1),nBins);
            [tempMedian2,~]=getDeciles(x5(tempIdx2),y5(tempIdx2),nBins);
            [tempMedian3,~]=getDeciles(x5(logical(tempIdx2.*tempIdx3)),y5(logical(tempIdx2.*tempIdx3)),nBins);
            [tempMedian4,~]=getDeciles(x5(logical(tempIdx2.*~tempIdx3)),y5(logical(tempIdx2.*~tempIdx3)),nBins);
            tempMedian5=(3*tempMedian3+tempMedian4)/4;
            
            subplot(4,4,8+l)
            hold on
            plot(1:length(tempMedian1),tempMedian1,'Color',blue)
            plot(1:length(tempMedian2),tempMedian2,'Color',orange)
            plot(1:length(tempMedian3),tempMedian3,'Color',grey)
            plot(1:length(tempMedian4),tempMedian4,':','Color','black')
            plot(1:length(tempMedian5),tempMedian5,'Color','black')
            xlabel(propertyLabels{xProperties(i)})
            ylabel(propertyLabels{yProperties(j)})
            ylim([yLim1(yProperties(j)) yLim2(yProperties(j))])
            title(['synonymous ' structureLabels{l}])
            xlim([1 nBins])
            axis square
            %legend({'observed','simulated','simulated Ts only','simulated Tv only','Ts/Tv=3'})
            
        end
        
                
        set(gcf,'PaperPositionMode','auto')
        print(['figures/structureSelection_' propertyLabels{xProperties(i)} '_' propertyLabels{yProperties(j)} '_fig_' num2str(figureCounter)],'-dsvg','-r0')
        print(['figures/structureSelection_' propertyLabels{xProperties(i)} '_' propertyLabels{yProperties(j)} '_fig_' num2str(figureCounter)],'-djpeg','-r0')
        figureCounter=figureCounter+1;

        close

        %rare vs common
        

        
        toc

    end
    
end


hrwgve


subplot(2,3,2)
vHasStruct=~cellfun(@isempty,secondary);
[v1,structArray]=structureTypes(secondary(vHasStruct));
bar(v1)
xticks(1:length(structureLabels))
xtickangle(45)
xticklabels(structureLabels)
title('in 1K genomes [all mis & syn]')
ylabel('frequency')
xticks(1:length(structureLabels))
xtickangle(45)
xticklabels(structureLabels)
ylim([0 0.5])


vMis=ismember(type,'missense_variant');
vSyn=ismember(type,'synonymous_variant');

[v1 misStructArray]=structureTypes(secondary(logical(vHasStruct.*vMis)));
[v2 synStructArray]=structureTypes(secondary(logical(vHasStruct.*vSyn)));
misSasaHasStruct=sasa(logical(vHasStruct.*vMis));
synSasaHasStruct=sasa(logical(vHasStruct.*vSyn));

misEssentialHasStruct=isEssential(logical(vHasStruct.*vMis));
synEssentialHasStruct=isEssential(logical(vHasStruct.*vSyn));

subplot(2,3,3)
bar([v1; v2]')
title('in 1K genomes')
legend({'mis','syn'})
ylabel('frequency')
xticks(1:length(structureLabels))
xtickangle(45)
xticklabels(structureLabels)
ylim([0 0.5])


for i=1:4
    v3(i)=sum(misSecondary==i);
    v4(i)=sum(synSecondary==i);
end
v3=v3./sum(v3);
v4=v4./sum(v4);

%only Ts and only Tv
for i=1:4
    v5(i)=sum(misSecondary(misTs==1)==i);
    v6(i)=sum(synSecondary(synTs==1)==i);
    
    v7(i)=sum(misSecondary(misTs==0)==i);
    v8(i)=sum(synSecondary(synTs==0)==i);
end
v5=v5./sum(v5);
v6=v6./sum(v6);

v7=v7./sum(v7);
v8=v8./sum(v8);

subplot(2,3,4)
bar([v1; v3; v5; v7]')
title('missense')
ylabel('frequency')
ylim([0 0.5])
xticks(1:length(structureLabels))
xtickangle(45)
xticklabels(structureLabels)
legend({'observed mis','simulated mis','simulated mis [only Ts]','simulated mis [only Tv]'})


subplot(2,3,5)
bar([v2; v4; v6; v8]')
title('synonymous')
ylabel('frequency')
ylim([0 0.5])
xticks(1:length(structureLabels))
xtickangle(45)
xticklabels(structureLabels)
legend({'observed syn','simulated syn','simulated syn [only Ts]','simulated syn [only Tv]'})




%actual vs prior sasa
clear toPlot
subplot(2,6,11)
hold on
toPlot{1}=sasa(vMis);
toPlot{2}=misSasa;
toPlot{3}=misSasa(misTs==1);
toPlot{4}=misSasa(misTs==0);
easyBox(toPlot)
title('missense')
ylabel('accessible surface area')
ylim([-50 300])
xticks(1:length(toPlot))
xtickangle(45)
xticklabels({'observed mis','simulated mis','simulated mis [only Ts]','simulated mis [only Tv]'})
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end
m=0;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [p h]=ranksum(toPlot{i},toPlot{j});
        plot([i j],[150+25*m 150+25*m],'-k')
        text((i+j)/2,160+25*m,num2str(p))
        m=m+1;
    end
end


subplot(2,6,12)
hold on
toPlot{1}=sasa(vSyn);
toPlot{2}=synSasa;
toPlot{3}=synSasa(synTs==1);
toPlot{4}=synSasa(synTs==0);
easyBox(toPlot)
title('synonymous')
ylabel('accessible surface area')
ylim([-50 300])
xticks(1:length(toPlot))
xtickangle(45)
xticklabels({'observed syn','simulated syn','simulated syn [only Ts]','simulated syn [only Tv]'})
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end
m=0;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [p h]=ranksum(toPlot{i},toPlot{j});
        plot([i j],[150+25*m 150+25*m],'-k')
        text((i+j)/2,160+25*m,num2str(p))
        m=m+1;
    end
end



set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;







%subset sasa (all vs observed) by structure type
clear toPlot
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    toPlot{m}=misSasaHasStruct(tempIdx1);
    tempLabels{m}=['observed mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=misSecondary==i;
    toPlot{m}=misSasa(tempIdx2);
    tempLabels{m}=['simulated mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==1;
    toPlot{m}=misSasa(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==0;
    toPlot{m}=misSasa(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1}))/4);
    
end
    


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([-50 200])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('accessible surface area')
title('missense variants')

subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('relative access. surface area (observed/simulated)')
title('missense variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')



clear toPlot
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    toPlot{m}=synSasaHasStruct(tempIdx1);
    tempLabels{m}=['observed syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=synSecondary==i;
    toPlot{m}=synSasa(tempIdx2);
    tempLabels{m}=['simulated syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==1;
    toPlot{m}=synSasa(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==0;
    toPlot{m}=synSasa(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1}))/4);
    
end
    


subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([-50 200])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('accessible surface area')
title('synonymous variants')

subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('relative access. surface area (observed/simulated)')
title('synonymous variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')



%SASA by af for observed
afHasStruct=af(vHasStruct);

%tempMedian=median(afHasStruct,'omitnan');
tempMedian=1e-3/2;    %rare [present only as 1 het]

afIdx=afHasStruct<tempMedian;

sasaHasStruct=sasa(vHasStruct);

clear toPlot
m=1;

for i=1:4
    tempIdx1=structArray==i;
    tempIdx2=afIdx;
    toPlot{m}=sasa(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed rare ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=structArray==i;
    tempIdx2=~afIdx;
    toPlot{m}=sasa(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed common ' structureLabels{i}];
    m=m+1;
    
end

% subplot(2,2,4)
% hold on
% easyBox(toPlot)
% ylim([-50 200])
% xticks(1:length(tempLabels))
% xtickangle(45)
% xticklabels(tempLabels)
% ylabel('accessible surface area')
% title('all observed variants')




set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;





%subset by missense vs synonymous

%SASA by af for observed
misAf=af(logical(vMis.*vHasStruct));
synAf=af(logical(vSyn.*vHasStruct));




%tempMedian=median(misAf,'omitnan');
tempMedian=1e-3/2;    %rare [present only as 1 het]

misAfIdx=misAf<=tempMedian;

clear tempLabels
clear toPlot
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    tempIdx2=misAfIdx;
    toPlot{m}=misSasaHasStruct(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed rare ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=misStructArray==i;
    tempIdx2=~misAfIdx;
    toPlot{m}=misSasaHasStruct(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed common ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
    
end



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([-50 200])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('accessible surface area')
title('missense variants by AF')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end

subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('relative access. surface area (rare/other)')
title('missense variants')
plot(xlim,[1 1],':k')





%tempMedian=median(synAf,'omitnan');
tempMedian=1e-3/2;    %rare [present only as 1 het]

synAfIdx=synAf<=tempMedian;

clear toPlot
clear tempLabels
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    tempIdx2=synAfIdx;
    toPlot{m}=synSasaHasStruct(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed rare ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=synStructArray==i;
    tempIdx2=~synAfIdx;
    toPlot{m}=synSasaHasStruct(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed common ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
end



subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([-50 200])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('accessible surface area')
title('synonymous variants by AF')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end


subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('relative access. surface area (rare/other)')
title('synonymous variants')
plot(xlim,[1 1],':k')




set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;




%analyze codon frequencies/transitions

[misInitialCodonFreq,misInitialnTE]=codonFrequency(misInitialCodon);
[misNewCodonFreq,misNewnTE]=codonFrequency(misNewCodon);

misDeltanTE=abs(misNewnTE-misInitialnTE);
misMeannTE=mean([misInitialnTE misNewnTE],2);

misInitialCodon1k=refCodon(logical(vMis.*vHasStruct));
misNewCodon1k=altCodon(logical(vMis.*vHasStruct));

[misInitialCodon1kFreq,misInitialCodon1knTE]=codonFrequency(misInitialCodon1k);
[misNewCodon1kFreq,misNewCodon1knTE]=codonFrequency(misNewCodon1k);

misDeltanTE1k=abs(misNewCodon1knTE-misInitialCodon1knTE);
misMeannTE1k=mean([misInitialCodon1knTE misNewCodon1knTE],2);

clear toPlot
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    toPlot{m}=misDeltanTE1k(tempIdx1);
    tempLabels{m}=['observed mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=misSecondary==i;
    toPlot{m}=misDeltanTE(tempIdx2);
    tempLabels{m}=['simulated mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==1;
    toPlot{m}=misDeltanTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==0;
    toPlot{m}=misDeltanTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1},'omitnan'))/4);
    
end
    


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('missense variants')


subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (observed/simulated)')
title('missense variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')




[synInitialCodonFreq,synInitialnTE]=codonFrequency(synInitialCodon);
[synNewCodonFreq,synNewnTE]=codonFrequency(synNewCodon);

synDeltanTE=abs(synNewnTE-synInitialnTE);
synMeannTE=mean([synInitialnTE synNewnTE],2);

synInitialCodon1k=refCodon(logical(vSyn.*vHasStruct));
synNewCodon1k=altCodon(logical(vSyn.*vHasStruct));

[synInitialCodon1kFreq,synInitialCodon1knTE]=codonFrequency(synInitialCodon1k);
[synNewCodon1kFreq,synNewCodon1knTE]=codonFrequency(synNewCodon1k);

synDeltanTE1k=abs(synNewCodon1knTE-synInitialCodon1knTE);
synMeannTE1k=mean([synInitialCodon1knTE synNewCodon1knTE],2);

clear toPlot
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    toPlot{m}=synDeltanTE1k(tempIdx1);
    tempLabels{m}=['observed syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=synSecondary==i;
    toPlot{m}=synDeltanTE(tempIdx2);
    tempLabels{m}=['simulated syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==1;
    toPlot{m}=synDeltanTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==0;
    toPlot{m}=synDeltanTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1},'omitnan'))/4);
    
end
    


subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('synonymous variants')


subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (observed/simulated)')
title('synonymous variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')


set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;




clear toPlot
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    toPlot{m}=misMeannTE1k(tempIdx1);
    tempLabels{m}=['observed mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=misSecondary==i;
    toPlot{m}=misMeannTE(tempIdx2);
    tempLabels{m}=['simulated mis. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==1;
    toPlot{m}=misMeannTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=misTs==0;
    toPlot{m}=misMeannTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated mis. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1},'omitnan'))/4);
    
end
    


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('mean nTE')
title('missense variants')


subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('mean nTE (observed/simulated)')
title('missense variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')



clear toPlot
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    toPlot{m}=synMeannTE(tempIdx1);
    tempLabels{m}=['observed syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx2=synSecondary==i;
    toPlot{m}=synMeannTE(tempIdx2);
    tempLabels{m}=['simulated syn. ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==1;
    toPlot{m}=synMeannTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Ts only] ' structureLabels{i}];
    m=m+1;
    
    tempIdx3=synTs==0;
    toPlot{m}=synMeannTE(logical(tempIdx2.*tempIdx3));
    tempLabels{m}=['simulated syn. [Tv only] ' structureLabels{i}];
    m=m+1;
    
    %let Ts/Tv=3
    oddsRatio(i)=median(toPlot{m-4},'omitnan')/...
        ((3*median(toPlot{m-2},'omitnan')+1*median(toPlot{m-1},'omitnan'))/4);
    
end
    


subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('mean nTE')
title('synonymous variants')


subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('mean nTE (observed/simulated)')
title('synonymous variants [let Ts/Tv=3]')
plot(xlim,[1 1],':k')




set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;






bases={'A','C','G','T'};

codons=cell(64,1);
m=1;
for i=1:length(bases)
    
    for j=1:length(bases)
        
        for k=1:length(bases)
            
            codons{m}=[bases{i} bases{j} bases{k}];
            m=m+1;
            
        end
        
    end
    
end


%which codon transitions are favored
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
hist3([refCodon,altCodon],[64 64],'CdataMode','auto')
view(2)
axis square
title('all 1k SNPs')

subplot(2,3,2)
hist3([misInitialCodon1k,misNewCodon1k],[64 64],'CdataMode','auto')
view(2)
axis square
title('1k missense')

subplot(2,3,3)
hist3([synInitialCodon1k,synNewCodon1k],[64 64],'CdataMode','auto')
view(2)
axis square
title('1k synonymous')





subplot(2,3,5)
hist3([[misInitialCodon(misTs==1);misInitialCodon(misTs==1);misInitialCodon(misTs==1);misInitialCodon(misTs==0)],...
    [misNewCodon(misTs==1);misNewCodon(misTs==1);misNewCodon(misTs==1);misNewCodon(misTs==0)]],[64 64],'CdataMode','auto')
view(2)
axis square
title('sim. missense [Ts/Tv=3]')

subplot(2,3,6)
hist3([[synInitialCodon(synTs==1);synInitialCodon(synTs==1);synInitialCodon(synTs==1);synInitialCodon(synTs==0)],...
    [synNewCodon(synTs==1);synNewCodon(synTs==1);synNewCodon(synTs==1);synNewCodon(synTs==0)]],[64 64],'CdataMode','auto')
view(2)
axis square
title('sim. synonymous [Ts/Tv=3]')



set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;






%tempMedian=median(misAf,'omitnan');
tempMedian=1e-3/2;    %rare [present only as 1 het]

misAfIdx=misAf<=tempMedian;

clear tempLabels
clear toPlot
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    tempIdx2=misAfIdx;
    toPlot{m}=misDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed rare ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=misStructArray==i;
    tempIdx2=~misAfIdx;
    toPlot{m}=misDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed common ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
    
end



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('missense variants by AF')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end

subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (rare/other)')
title('missense variants')
plot(xlim,[1 1],':k')





%tempMedian=median(synAf,'omitnan');
tempMedian=1e-3/2;    %rare [present only as 1 het]

synAfIdx=synAf<=tempMedian;

clear toPlot
clear tempLabels
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    tempIdx2=synAfIdx;
    toPlot{m}=synDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed rare ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=synStructArray==i;
    tempIdx2=~synAfIdx;
    toPlot{m}=synDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed common ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
end



subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('synonymous variants by AF')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end


subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (rare/other)')
title('synonymous variants')
plot(xlim,[1 1],':k')




set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;



%delta nTE by sasa (above or below median) for observed
tempMedian=median(misSasaHasStruct);


clear toPlot
clear tempLabels
m=1;

for i=1:4
    tempIdx1=misStructArray==i;
    tempIdx2=misSasaHasStruct>tempMedian;
    toPlot{m}=misDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed exposed ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=misStructArray==i;
    tempIdx2=misSasaHasStruct<tempMedian;
    toPlot{m}=misDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed buried ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
end


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('missense variants by acc. surface area')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end



subplot(2,4,5)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (exposed/buried)')
title('missense variants')
plot(xlim,[1 1],':k')




%delta nTE by sasa (above or below median) for observed
tempMedian=median(synSasaHasStruct);


clear toPlot
clear tempLabels
m=1;

for i=1:4
    tempIdx1=synStructArray==i;
    tempIdx2=synSasaHasStruct>tempMedian;
    toPlot{m}=synDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed exposed ' structureLabels{i}];
    m=m+1;
    
    tempIdx1=synStructArray==i;
    tempIdx2=synSasaHasStruct<tempMedian;
    toPlot{m}=synDeltanTE1k(logical(tempIdx1.*tempIdx2'));
    tempLabels{m}=['observed buried ' structureLabels{i}];
    m=m+1;
    
    oddsRatio(i)=median(toPlot{m-2},'omitnan')/median(toPlot{m-1},'omitnan');
    
end


subplot(2,2,2)
hold on
easyBox(toPlot)
ylim([0 0.5])
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
ylabel('|\DeltanTE|')
title('synonymous variants by acc. surface area')
for i=1:(length(tempLabels)/2)
    [p h]=ranksum(toPlot{2*(i-1)+1},toPlot{2*(i-1)+2});
    plot([2*(i-1)+1 2*(i-1)+2],[-35 -35],'-k')
    text(2*(i-1)+1.5,-25,num2str(p))
end
for i=1:length(toPlot)
    text(i,-40,num2str(sum(~isnan(toPlot{i}))))
end



subplot(2,4,6)
hold on
bar(oddsRatio)
ylim([0.5 2])
set(gca,'YScale','log')
xticks(1:4)
xtickangle(45)
xticklabels(structureLabels)
ylabel('|\DeltanTE| (exposed/buried)')
title('synonymous variants')
plot(xlim,[1 1],':k')





set(gcf,'PaperPositionMode','auto')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-dsvg','-r0')
print(['figures/structureSelection_fig_' num2str(figureCounter)],'-djpeg','-r0')
figureCounter=figureCounter+1;





%codon spectrum
figure('units','normalized','outerposition',[0 0 1 1])
hold on
v1=misInitialCodon1kFreq./sum(misInitialCodon1kFreq);
v2=misNewCodon1kFreq./sum(misNewCodon1kFreq);
v3=misInitialCodonFreq./sum(misInitialCodonFreq);
v4=misNewCodonFreq./sum(misNewCodonFreq);
bar([v1 v2 v3 v4]) 
xticks(1:length(codons))
xtickangle(45)
xticklabels(codons)
title('missense')


figure('units','normalized','outerposition',[0 0 1 1])
hold on
v1=synInitialCodon1kFreq./sum(synInitialCodon1kFreq);
v2=synNewCodon1kFreq./sum(synNewCodon1kFreq);
v3=synInitialCodonFreq./sum(synInitialCodonFreq);
v4=synNewCodonFreq./sum(synNewCodonFreq);
bar([v1 v2 v3 v4]) 
xticks(1:length(codons))
xtickangle(45)
xticklabels(codons)
title('synonymous')


