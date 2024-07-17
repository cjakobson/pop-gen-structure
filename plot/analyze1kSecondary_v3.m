%condensed version of structure analysis


clear

tic

figureCounter=1; %figure counter

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)

filebase='/Users/cjakobson/Dropbox/JaroszLab/';


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


% structureLabels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
%     'bend','turn','unstr.'};
% structureLabels={'helices','sheet','ext. str.','turns','unstruct.'};
structureLabels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

doSim=1;





%test essential vs not
essentialList=readtable('Giaever_G_et_al_2002_phenotypes_filtered_by_inv.txt');
essentialList=table2array(essentialList(:,2));
%check essentiality
isEssential=ismember(gene,essentialList);

%idxToUse=isEssential;


idxToUse=logical(ones(length(gene),1));




load('secondary1k.mat','secondary')
load('sasa1k.mat','sasa')
load('neighbors1k.mat','neighbors')



load('simulationStructureData.mat')

%clear out some unnecessary variables
clear chr dnaEncoded dnaSequence dsspTable inputGenome %gene ...
clear mutationTable outputTable pdbDictionary pdbId

%organize data to plot systematically
propertyLabels={'accessible surface area','neighbors'};

yProperties=1:length(propertyLabels);
xProperties=yProperties;

yLim1=[0 0];
%yLim2=[150 2 180 180 1 50 50 1 0.5 0.5 0.5 25 50];
yLim2=[250 40];

type=type(idxToUse);
vMis1k=ismember(type,'missense_variant');

clear type

secondary=secondary(idxToUse);
vHasStruct=~cellfun(@isempty,secondary);

[~,structureMis1k]=structureTypes(secondary(logical(vMis1k.*vHasStruct)));  %for categorizing

clear secondary

sasa=sasa(idxToUse);
neighbors=neighbors(idxToUse);

properties1k{1,1}=sasa(logical(vMis1k.*vHasStruct));
properties1k{1,2}=neighbors(logical(vMis1k.*vHasStruct));


%allele frequencies
af=af(idxToUse);
misAf=af(logical(vMis1k.*vHasStruct));


clear af sasa phi psi posInOrf posInRun lengthOfRun neighbors pDist %altCodon refCodon


%structToUse=[1:2 5:8];
structToUse=1:length(structureLabels);

structureMisSim=misSecondary;   %for categorizing

propertiesSim{1,1}=misSasa;
propertiesSim{1,2}=misNeighbors;

clear misSasa misPhi misPsi misPosInOrf misPosInRun misLengthOfRun misNeighbors misInitialCodon misNewCodon 
clear synSasa synPhi synPsi synPosInOrf synPosInRun synLengthOfRun synNeighbors %synInitialCodon synNewCodon 



%all segregating and trans missense pQTNs to compare
variantInfo=readtable([filebase '211028_SegregantProteomicsData_V1/variantInfoStructure.csv']);

allPqtl=readtable([filebase '211028_SegregantProteomicsData_V1/odCovariate/linearPqtlOd_FDR_0.1.csv']);

%remove OD600 term
allPqtl(allPqtl.bPos==1,:)=[];


propertiesAllSegregating{1}=variantInfo.variantSasa(ismember(variantInfo.variantType,'missense_variant'));
propertiesAllSegregating{2}=variantInfo.variantNeighbors(ismember(variantInfo.variantType,'missense_variant'));




qtnIdx=allPqtl.isQtn==1;
misIdx=ismember(allPqtl.variantType,'missense_variant');

pqtnIndex=unique(allPqtl.index(logical(qtnIdx.*misIdx)));

propertiesPqtn{1}=variantInfo.variantSasa(pqtnIndex);
propertiesPqtn{2}=variantInfo.variantNeighbors(pqtnIndex);




misIdx=ismember(variantInfo.variantType,'missense_variant');
allOtherIndex=variantInfo.index(misIdx);
allOtherIndex(ismember(allOtherIndex,pqtnIndex))=[];

propertiesAllOther{1}=variantInfo.variantSasa(allOtherIndex);
propertiesAllOther{2}=variantInfo.variantNeighbors(allOtherIndex);





toc



afThresh=1e-3;


for i=1:length(propertyLabels)

    tic
    
    %gross distribution (before considering stucture)
    figure('units','normalized','outerposition',[0 0 1 1])

    
    %mis then syn; merge to make all
    v1=properties1k{1,i};
    
    v4=propertiesSim{1,i};
    
    
    %compare all segregating to all possible
    clear toPlot
    
    %missense
    toPlot{1}=v4;
    toPlot{2}=propertiesAllSegregating{i};
    
    %histograms
    subplot(2,8,1)
    hold on
    easyBox(toPlot)
    ylim([0 yLim2(i)])
    xticklabels({'all poss.','all segr.'})
    
    
    
    
    
    %also pQTNs to all other segregating
    clear toPlot
    
    %missense
    toPlot{1}=v4;
    toPlot{2}=propertiesAllOther{i};
    toPlot{3}=propertiesPqtn{i};
    
    %histograms
    subplot(2,8,2)
    hold on
    easyBox(toPlot)
    ylim([0 yLim2(i)])
    xticklabels({'all poss.','all other','pQTNs'})
    
    
    
    
    clear toPlot
    
    %missense
    toPlot{1}=v1;
    toPlot{2}=v4;
    
    %histograms
    subplot(2,4,3)
    hold on
    histogram(toPlot{1},0:(ceil(yLim2(i)/20)):(yLim2(i)),'Normalization','Probability')
    histogram(toPlot{2},0:(ceil(yLim2(i)/20)):(yLim2(i)),'Normalization','Probability')
    xlim([0 yLim2(i)])
    %xlim([0 500])
    axis square
    title('missense')
    legend({'1K','all poss.'})
    [p h]=ranksum(toPlot{1},toPlot{2});
    text(40,0.08,num2str(p))
    
    
    
    clear toPlot tempLabels pVal
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        
        %sim
        tempIdx2=misSecondary==j;
        
        %toPlot(m)=median(v1(tempIdx1),'omitnan')/median(v4(tempIdx2),'omitnan');
        toPlot(m)=mean(v1(tempIdx1),'omitnan')/mean(v4(tempIdx2),'omitnan');
        tempLabels{m}=['obs./sim. ' structureLabels{structToUse(j)}];
        [pVal(m) h]=ranksum(v1(tempIdx1),v4(tempIdx2));
        m=m+1;
        
        
    end
    
    subplot(2,6,7)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel([propertyLabels{i} ' 1K/all poss.'])
    title('missense')
    xlim([0.5 length(tempLabels)+0.5])
    ylim([0.8 1.3])
    %ylim([1/10 10])
    set(gca,'YScale','log')
    for j=1:length(pVal)
        text(j,1.3,num2str(pVal(j)))
    end
    
    
    %rare vs common odds ratio plots?
    clear toPlot
    %figure('units','normalized','outerposition',[0 0 1 1])

    %mis then syn; merge to make all
    v1=properties1k{1,i};
    
    misIdx=misAf<afThresh;
    
    toPlot{1}=v1(misIdx);
    toPlot{2}=v1(~misIdx);
    
    
    %also histogram
    subplot(2,4,4)
    hold on
    histogram(toPlot{1},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    histogram(toPlot{2},0:(ceil(yLim2(i)/25)):(2*yLim2(i)),'Normalization','Probability')
    xlim([0 yLim2(i)])
    axis square
    title('missense')
    legend({'rare','common'})
    [p h]=ranksum(toPlot{1},toPlot{2});
    text(40,0.08,num2str(p))
    
    
    
    
    %subset by structure
    clear toPlot
    clear tempLabels
    
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        tempIdx2=misIdx;
        
        
%         toPlot(m)=median(v1(logical(tempIdx1.*~tempIdx2)),'omitnan')/...
%             median(v1(logical(tempIdx1.*tempIdx2)),'omitnan');
        toPlot(m)=mean(v1(logical(tempIdx1.*~tempIdx2)),'omitnan')/...
            mean(v1(logical(tempIdx1.*tempIdx2)),'omitnan');
        tempLabels{m}=['rare/common ' structureLabels{structToUse(j)}];
        [pVal(m) h]=ranksum(v1(logical(tempIdx1.*~tempIdx2)),...
            v1(logical(tempIdx1.*tempIdx2)));
        m=m+1;
        
    end
    
    subplot(2,6,8)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel([propertyLabels{i} ' common/rare'])
    xlim([0.5 length(tempLabels)+0.5])
    %ylim([2/3 3/2])
    ylim([0.85 1.2])
    set(gca,'YScale','log')
    title('missense 1K')
    %pVals
    for j=1:length(pVal)
        text(j,1.15,num2str(pVal(j)))
    end
    
    
    
    
    %subset by essential vs not within 1K
    tempIdx3=isEssential(logical(vMis1k.*vHasStruct));
    
    clear toPlot
    clear tempLabels
    
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        tempIdx2=misIdx;
        
        
%         toPlot(m)=median(v1(logical(tempIdx1.*~tempIdx2)),'omitnan')/...
%             median(v1(logical(tempIdx1.*tempIdx2)),'omitnan');
        toPlot(m)=mean(v1(logical(tempIdx1.*~tempIdx2.*tempIdx3)),'omitnan')/...
            mean(v1(logical(tempIdx1.*tempIdx2.*tempIdx3)),'omitnan');
        tempLabels{m}=['rare/common ' structureLabels{structToUse(j)}];
        [pVal(m) h]=ranksum(v1(logical(tempIdx1.*~tempIdx2.*tempIdx3)),...
            v1(logical(tempIdx1.*tempIdx2.*tempIdx3)));
        m=m+1;
        
    end
    
    subplot(2,6,9)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel([propertyLabels{i} ' common/rare'])
    xlim([0.5 length(tempLabels)+0.5])
    %ylim([2/3 3/2])
    ylim([0.85 1.2])
    set(gca,'YScale','log')
    title('missense 1K essential genes')
    %pVals
    for j=1:length(pVal)
        text(j,1.15,num2str(pVal(j)))
    end
    
    
    
    clear toPlot
    clear tempLabels
    
    m=1;
    for j=1:length(structureLabels(structToUse))
        
        %observed
        tempIdx1=structureMis1k'==j;
        tempIdx2=misIdx;
        
        
%         toPlot(m)=median(v1(logical(tempIdx1.*~tempIdx2)),'omitnan')/...
%             median(v1(logical(tempIdx1.*tempIdx2)),'omitnan');
        toPlot(m)=mean(v1(logical(tempIdx1.*~tempIdx2.*~tempIdx3)),'omitnan')/...
            mean(v1(logical(tempIdx1.*tempIdx2.*~tempIdx3)),'omitnan');
        tempLabels{m}=['rare/common ' structureLabels{structToUse(j)}];
        [pVal(m) h]=ranksum(v1(logical(tempIdx1.*~tempIdx2.*~tempIdx3)),...
            v1(logical(tempIdx1.*tempIdx2.*~tempIdx3)));
        m=m+1;
        
    end
    
    subplot(2,6,10)
    hold on
    bar(toPlot,'BaseValue',1)
    xticks(1:length(tempLabels))
    xtickangle(45)
    xticklabels(tempLabels)
    ylabel([propertyLabels{i} ' common/rare'])
    xlim([0.5 length(tempLabels)+0.5])
    %ylim([2/3 3/2])
    ylim([0.85 1.2])
    set(gca,'YScale','log')
    title('missense 1K non-essential genes')
    %pVals
    for j=1:length(pVal)
        text(j,1.15,num2str(pVal(j)))
    end
    
    
    
    
    set(gcf,'PaperPositionMode','auto')
    print(['figures/structureSelection_' propertyLabels{i}  '_fig_' num2str(figureCounter)],'-dsvg','-r0')
    print(['figures/structureSelection_' propertyLabels{i}  '_fig_' num2str(figureCounter)],'-djpeg','-r0')
    figureCounter=figureCounter+1;
    
    
    %close
    
    
    toc
    
    
    
end


