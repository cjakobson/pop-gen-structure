%analyze primordium sequencing of FPR1


clear

tic

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

%filebase='/Users/cjakobson/';
filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/lpop-gen-structure-output/'];


addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])
addpath([code_directory 'data-prep'])

to_read=dir([dependency_directory 'FPR1_results/K2P4DZ_per_base_data']);
to_read=to_read(3:end);



for i=1:length({to_read.name})
    
    temp_str=strsplit(to_read(i).name,'_');
    v_to_sort(i)=str2num(temp_str{2}(1:(end-4)));
    
end

[v_sort,sort_idx2]=sort(v_to_sort,'ascend');



base_array={'A','C','T','G'};

read_thresh=40;%100;
fraction_thresh=0.33;

missense_position=cell(length(to_read),1);
ref_base=missense_position;
alt_base=missense_position;
missense_encoded=missense_position;
indel_position=missense_position;
indel_base=missense_position;

all_fraction=[];
for i=1:length(to_read)
    
    input_table=readtable([dependency_directory 'FPR1_results/K2P4DZ_per_base_data/' to_read(sort_idx2(i)).name]);
    
    input_table=input_table(33:end,:);
    
    v_mismatch=table2array(input_table(:,5));
    v_counts=table2array(input_table(:,3));

    v_fraction=v_mismatch./v_counts;

    all_fraction=[all_fraction;v_fraction];
    
    %temp_idx=find(v_mismatch>read_thresh);
    temp_idx=find(v_fraction*12>=fraction_thresh);
    
    
    for j=1:length(temp_idx)
        
        [v_temp,sort_idx]=sort(table2array(input_table(temp_idx(j),8:11)),'descend');
        
        missense_position{i}(j)=temp_idx(j);
        ref_base{i}{j}=base_array{sort_idx(1)};
        alt_base{i}{j}=base_array{sort_idx(2)};
        missense_fraction{i}(j)=v_fraction(temp_idx(j));
        
    end
    
    v_indel=max(table2array(input_table(:,6:7)),[],2);
    
    temp_idx=find(v_indel>read_thresh);
    
    
    for j=1:length(temp_idx)
        
        [v_temp,sort_idx]=sort(table2array(input_table(temp_idx(j),8:11)),'descend');
        
        indel_position{i}(j)=temp_idx(j);
        indel_base{i}{j}=base_array{sort_idx(1)};
        
    end
    
end


%only use rap no rad pools
pools_to_use=[1:3 7:9 13:15 19:21];









%need to account for direction and offset relative to FPR1 start codon
v_direction={'fwd','rev','fwd','fwd','rev','rev',...
    'rev','fwd','fwd','fwd','fwd','fwd',...
    'fwd','fwd','rev','fwd','rev','fwd',...
    'rev','rev','rev','rev','fwd','rev'};

fpr1_length=345;

missense_pos_fpr1=cell(length(to_read),1);
indel_pos_fpr1=cell(length(to_read),1);
for i=1:length(missense_position)
    
    if strcmp(v_direction{i},'fwd')
        
        v_offset1=216;%215;
        
        for j=1:length(missense_position{i})
            
            temp_pos=missense_position{i}(j);
            
            missense_pos_fpr1{i}(j)=temp_pos-v_offset1;
            
        end
        
        for j=1:length(indel_position{i})
            
            temp_pos=indel_position{i}(j);
            
            indel_pos_fpr1{i}(j)=temp_pos-v_offset1;
            
        end
        
    end
    
    if strcmp(v_direction{i},'rev')
        
        v_offset2=235;%234;
        
        for j=1:length(missense_position{i})
            
            temp_pos=missense_position{i}(j);
            
            temp_pos2=temp_pos-v_offset2;
            
            missense_pos_fpr1{i}(j)=fpr1_length-temp_pos2+1;
            
            %rev complement of bases
            ref_base{i}{j}=seqrcomplement(ref_base{i}{j});
            alt_base{i}{j}=seqrcomplement(alt_base{i}{j});
            
        end
        
        for j=1:length(indel_position{i})
            
            temp_pos=indel_position{i}(j);
            
            temp_pos2=temp_pos-v_offset2;
            
            indel_pos_fpr1{i}(j)=fpr1_length-temp_pos2+1;
            
        end
        
    end
    
end


fpr1_seq=fastaread([dependency_directory 'S288C_YNL135C_FPR1_coding.fsa']);

for i=1:length(missense_encoded)
    
    for j=1:length(missense_pos_fpr1{i})
        
        temp_residue=ceil(missense_pos_fpr1{i}(j)/3);
        
        temp_range=(3*temp_residue-2):(3*temp_residue);
        
        if (temp_range(1)>0)&&(temp_range(2)<length(fpr1_seq.Sequence))

            temp_ref_codon=fpr1_seq.Sequence(temp_range);
            
            temp_idx=ismember(temp_range,missense_pos_fpr1{i}(j));
            
            temp_alt_codon=temp_ref_codon;
            temp_alt_codon(temp_idx)=alt_base{i}{j};
            
            missense_encoded{i}{j}=[nt2aa(temp_ref_codon,'AlternativeStartCodons','false')...
                num2str(temp_residue) nt2aa(temp_alt_codon,'AlternativeStartCodons','false')];

        else

            missense_encoded{i}{j}='*';

        end
        
    end
    
end

missense_position=missense_position(pools_to_use);
missense_encoded=missense_encoded(pools_to_use);
missense_fraction=missense_fraction(pools_to_use);

indel_position=indel_position(pools_to_use);
indel_pos_fpr1=indel_pos_fpr1(pools_to_use);

v_direction=v_direction(pools_to_use);

%output table by pool
m=1;
for i=1:length(pools_to_use)

    pool_names{i}=['pool_' num2str(i)];

    if strcmp(v_direction{i},'fwd')
    
        [v_sort,sort_idx]=sort(missense_position{i},'ascend');

    elseif strcmp(v_direction{i},'rev')
    
        [v_sort,sort_idx]=sort(missense_position{i},'descend');

    end

    v_temp=[];
    for j=1:length(missense_encoded{i})

        v_temp=[v_temp missense_encoded{i}{sort_idx(j)} '_'...
            num2str(12*missense_fraction{i}(sort_idx(j))) ' '];

    end

    pool_mutations{i}=v_temp;

    %filter stops and count missense
    for j=1:length(missense_encoded{i})

        if ~strcmp(missense_encoded{i}{j}(end),'*')

            all_missense_identified{m}=missense_encoded{i}{j};
            m=m+1;

        end

    end

end

to_output=table(pool_names',pool_mutations',...
    'VariableNames',{'pool','mutations'});
writetable(to_output,[dependency_directory 'FPR1_primordium_SNPs.txt'])

all_missense_identified=sort(all_missense_identified);
unique_missense_identified=unique(all_missense_identified);


sum(cellfun(@length,missense_encoded))



bevdc







%calculate various structure statistics from alphafold calcs


load('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/1002data/1002dataAnn.mat')

temp_idx=cellfun(@isempty,gene);

%only need gene, position, variant type, and AF for now
gene(temp_idx)=[];
proteinEncoded(temp_idx)=[];
type(temp_idx)=[];
af(temp_idx)=[];


orfToUse='YNL135C';

load(['/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'neighborTable.mat'])
neighborTable=outputTable;

load(['/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'dsspTable.mat'])
dsspTable=outputTable;

load(['/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'mutationTable.mat'])


simulatedResidue=outputTable.residueNumber(outputTable.isMis==1);

geneIdx=ismember(gene,orfToUse);
typeIdx=ismember(type,'missense_variant');

observedMutations=proteinEncoded(logical(geneIdx.*typeIdx));

clear chr dnaEncoded gene pos proteinEncoded type

for i=1:length(observedMutations)
    
    observedResidue(i)=str2num(observedMutations{i}(6:(end-3)));
    
end


for i=1:max(simulatedResidue)
    
    tempSim=sum(simulatedResidue==i);
    tempObs=sum(observedResidue==i);
    
    vRatio(i)=tempObs/tempSim;
    
end





%no rad
allNoRad=[];
for i=[1:3 7:9 13:15 19:21]%1:2
    
    for j=1:length(missense_encoded{i})
    
        if ~strcmp(missense_encoded{i}{j}(end),'*')
    
            allNoRad=[allNoRad missense_pos_fpr1{i}(j)];
            
        end
        
    end
    
end

%allNoRad=unique(allNoRad);

%convert to residue
for i=1:length(allNoRad)
    
    allNoRad(i)=ceil(allNoRad(i)/3);
    
end


allRad=[];
for i=[4:6 10:12 16:18 22:24]%3:4
     
    for j=1:length(missense_encoded{i})
    
        if ~strcmp(missense_encoded{i}{j}(end),'*')
    
            allRad=[allRad missense_pos_fpr1{i}(j)];
            
        end
        
    end
    
end


%allRad=unique(allRad);

%convert to residue
for i=1:length(allRad)
    
    allRad(i)=ceil(allRad(i)/3);
    
end





figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,8,1)
hold on

clear toPlot
toPlot{1}=dsspTable.sasa(allNoRad);
toPlot{2}=dsspTable.sasa;

easyBox(toPlot)
ylim([0 200])
title('SASA')
tempLabels={'resistant no rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [h p]=ttest2(toPlot{i},toPlot{j});
        text((i+j)/2,100+25*m,num2str(p))
        m=m+1;
    end
end


subplot(2,8,2)
hold on

clear toPlot
toPlot{1}=neighborTable.neighbors(allNoRad);
toPlot{2}=neighborTable.neighbors;

easyBox(toPlot)
ylim([0 30])
title('neighbors')
tempLabels={'resistant no rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [h p]=ttest2(toPlot{i},toPlot{j});
        text((i+j)/2,25+1*m,num2str(p))
        m=m+1;
    end
end





subplot(2,8,9)
hold on

clear toPlot
toPlot{1}=dsspTable.sasa(allNoRad);
toPlot{2}=dsspTable.sasa(allRad);
toPlot{3}=dsspTable.sasa;

easyBox(toPlot)
ylim([0 200])
title('SASA')
tempLabels={'resistant no rad','resistant rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [h p]=ttest2(toPlot{i},toPlot{j});
        text((i+j)/2,100+25*m,num2str(p))
        m=m+1;
    end
end




subplot(2,8,10)
hold on

clear toPlot
toPlot{1}=neighborTable.neighbors(allNoRad);
toPlot{2}=neighborTable.neighbors(allRad);
toPlot{3}=neighborTable.neighbors;

easyBox(toPlot)
ylim([0 30])
title('neighbors')
tempLabels={'resistant no rad','resistant rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(toPlot)
    for j=(i+1):length(toPlot)
        [h p]=ttest2(toPlot{i},toPlot{j});
        text((i+j)/2,25+1*m,num2str(p))
        m=m+1;
    end
end






set(gcf,'PaperPositionMode','auto')
print(['fpr1_mutations_' num2str(figureCounter)],'-dsvg','-r0')
print(['fpr1_mutations_' num2str(figureCounter)],'-djpeg','-r300')
figureCounter=figureCounter+1;






