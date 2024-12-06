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

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

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
fraction_thresh_missense=0.33;
fraction_thresh_indel=2*fraction_thresh_missense;%0.8;

n_clones_in_pool=[12 12 12 12 12 9,...
    12 12 12 12 12 9,...
    12 12 8 12 12 12,...
    12 12 7 12 12 12];

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
    temp_idx=find(v_fraction*n_clones_in_pool(i)>=fraction_thresh_missense);
    
    
    for j=1:length(temp_idx)
        
        [v_temp,sort_idx]=sort(table2array(input_table(temp_idx(j),8:11)),'descend');
        
        missense_position{i}(j)=temp_idx(j);
        ref_base{i}{j}=base_array{sort_idx(1)};
        alt_base{i}{j}=base_array{sort_idx(2)};
        missense_fraction{i}(j)=v_fraction(temp_idx(j));
        
    end
    
    v_indel=max(table2array(input_table(:,6:7)),[],2);
    
    
    v_fraction=v_indel./v_counts;

    
    
    %temp_idx=find(v_mismatch>read_thresh);
    temp_idx=find(v_fraction*n_clones_in_pool(i)>=fraction_thresh_indel);

    
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
        
        v_offset1=215;%216;
        
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
        
        v_offset2=234;%235;
        
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

n_clones_in_pool=n_clones_in_pool(pools_to_use);

%filter for indels in gene
length_thresh=fpr1_length;
for i=1:length(indel_pos_fpr1)
    
    temp_idx=logical((indel_pos_fpr1{i}>0).*(indel_pos_fpr1{i}<=length_thresh));
    
    indel_position{i}=indel_position{i}(temp_idx);
    indel_pos_fpr1{i}=indel_pos_fpr1{i}(temp_idx);
    
    %merge neighboring (same indel)
    to_clear=zeros(1,length(indel_pos_fpr1{i}));
    for j=2:length(indel_pos_fpr1{i})
        
        if abs(indel_pos_fpr1{i}(j)-indel_pos_fpr1{i}(j-1))==1
            
            to_clear(j)=1;
            
        end
        
    end
    
    indel_position{i}(logical(to_clear))=[];
    indel_pos_fpr1{i}(logical(to_clear))=[];

    
end

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
            num2str(n_clones_in_pool(i)*missense_fraction{i}(sort_idx(j))) ' '];

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


%make output tables for missense, LoF, and indels
m=1;
for i=1:length(missense_encoded)
    
    temp_array=missense_encoded{i};
    
    for j=1:length(temp_array)
        
        temp_str=temp_array{j};
        
        snp_encoded{m}=temp_str;
        
        snp_ref{m}=temp_str(1);
        snp_alt{m}=temp_str(end);
        
        snp_pos(m)=str2num(temp_str(2:(end-1)));
        
        snp_frequency(m)=n_clones_in_pool(i)*missense_fraction{i}(j);
        
        snp_pool(m)=i;
        
        m=m+1;
        
    end
    
end

%exclude start codon and premature stops
temp_idx1=snp_pos>1;
temp_idx2=~ismember(snp_alt,'*');

to_use=logical(temp_idx1.*temp_idx2);

missense_pos_to_output=snp_pos(to_use);
missense_encoded_to_output=snp_encoded(to_use);
missense_pool_to_output=snp_pool(to_use);
missense_frequency_to_output=snp_frequency(to_use);

clone_thresh=1.5;
missense_inferred_clones=missense_frequency_to_output;
missense_inferred_clones(missense_inferred_clones>clone_thresh)=2;
missense_inferred_clones(missense_inferred_clones<=clone_thresh)=1;

[~,sort_idx]=sort(missense_pos_to_output,'ascend');

to_output=table(missense_pos_to_output(sort_idx)',...
    missense_encoded_to_output(sort_idx)',...
    missense_pool_to_output(sort_idx)',missense_frequency_to_output(sort_idx)',...
    missense_inferred_clones(sort_idx)');
writetable(to_output,[dependency_directory 'FPR1_primordium_missense.csv'])


lof_pos_to_output=snp_pos(~to_use);
lof_encoded_to_output=snp_encoded(~to_use);
lof_pool_to_output=snp_pool(~to_use);
lof_frequency_to_output=snp_frequency(~to_use);

lof_inferred_clones=lof_frequency_to_output;
lof_inferred_clones(lof_inferred_clones>clone_thresh)=2;
lof_inferred_clones(lof_inferred_clones<=clone_thresh)=1;

[~,sort_idx]=sort(lof_pos_to_output,'ascend');

to_output=table(lof_pos_to_output(sort_idx)',...
    lof_encoded_to_output(sort_idx)',...
    lof_pool_to_output(sort_idx)',lof_frequency_to_output(sort_idx)',...
    lof_inferred_clones(sort_idx)');
writetable(to_output,[dependency_directory 'FPR1_primordium_lof.csv'])


hbrva



%calculate various structure statistics from alphafold calcs


load([filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/1002data/1002dataAnn.mat'])

temp_idx=cellfun(@isempty,gene);

%only need gene, position, variant type, and AF for now
gene(temp_idx)=[];
proteinEncoded(temp_idx)=[];
type(temp_idx)=[];
af(temp_idx)=[];


orfToUse='YNL135C';

load([filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'neighborTable.mat'])
neighbor_table=outputTable;

load([filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'dsspTable.mat'])
dssp_table=outputTable;

load([filebase 'Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/mutationSimulations/mutationOutput/' orfToUse 'mutationTable.mat'])


simulated_sesidue=outputTable.residueNumber(outputTable.isMis==1);

gene_idx=ismember(gene,orfToUse);
type_idx=ismember(type,'missense_variant');

observed_mutations=proteinEncoded(logical(gene_idx.*type_idx));

clear chr dnaEncoded gene pos proteinEncoded type

for i=1:length(observed_mutations)
    
    observed_residue(i)=str2num(observed_mutations{i}(6:(end-3)));
    
end


for i=1:max(simulated_sesidue)
    
    temp_sim=sum(simulated_sesidue==i);
    temp_obs=sum(observed_residue==i);
    
    v_ratio(i)=temp_obs/temp_sim;
    
end


%merge missense variant positions
all_no_rad=[];
for i=1:length(missense_encoded)
    
    v_temp=missense_encoded{i};
    
    for j=1:length(v_temp)
        
        if ~strcmp(v_temp{j}(end),'*')
            
            all_no_rad=[all_no_rad str2num(v_temp{j}(2:(end-1)))];
            
        end
        
    end
    
end



figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,8,1)
hold on

clear toPlot
to_plot{1}=dssp_table.sasa(all_no_rad);
to_plot{2}=dssp_table.sasa;

easy_box(to_plot)
ylim([0 200])
title('SASA')
tempLabels={'resistant no rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(to_plot)
    for j=(i+1):length(to_plot)
        [h p]=ttest2(to_plot{i},to_plot{j});
        text((i+j)/2,100+25*m,num2str(p))
        m=m+1;
    end
end


subplot(2,8,2)
hold on

clear toPlot
to_plot{1}=neighbor_table.neighbors(all_no_rad);
to_plot{2}=neighbor_table.neighbors;

easy_box(to_plot)
ylim([0 30])
title('neighbors')
tempLabels={'resistant no rad','all'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
m=1;
for i=1:length(to_plot)
    for j=(i+1):length(to_plot)
        [h p]=ttest2(to_plot{i},to_plot{j});
        text((i+j)/2,25+1*m,num2str(p))
        m=m+1;
    end
end







