
clear

tic

%filebase='/Users/cjakobson/';
filebase='/Users/christopherjakobson/';

figure_counter=1;

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-output/'];


addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

load([dependency_directory 'asa_data.mat'])
load([dependency_directory 'neighbor_data.mat'])
load([dependency_directory 'residue_data.mat'])
load([dependency_directory 'secondary_structure_data.mat'])
load([dependency_directory 'superfam_domain_boundaries.mat'])

%make compatible with ismember
residue_mat(cellfun(@isempty,residue_mat))={'NA'};
secondary_mat(cellfun(@isempty,secondary_mat))={'NA'};

load([dependency_directory 'simulated_mutation_positions.mat'])

load([dependency_directory '1K_mutation_positions.mat'])
load([dependency_directory '1K_mutation_af.mat'])

%change to maf
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);


structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


hap_table=readtable([dependency_directory 'pbio.2005130.s003.xlsx'],'Sheet','A');

essential_table=readtable([dependency_directory 'inviable_annotations_filtered_by_giaever.txt']);


%general histograms/correlations for structure stats (just for all
%residues)
figure('units','normalized','outerposition',[0 0 1 1])
v1=reshape(asa_mat,1,[]);
v1=v1(~isnan(v1));
v2=reshape(neighbor_mat,1,[]);
v2=v2(~isnan(v2));

subplot(2,4,1)
hold on
histogram(v1,0:10:250)%,'Normalization','probability')
axis square
legend({'all residues'})
title('ASA')
ylabel('number of residues')
set(gca,'YScale','log')


subplot(2,4,2)
hold on
histogram(v2,0:2:50)%,'Normalization','probability')
axis square
legend({'all residues'})
title('neighbors')
ylabel('number of residues')
set(gca,'YScale','log')



subplot(2,4,3)
%scatter(v1,v2,10,'k','filled')
histogram2(v1,v2,0:5:250,4:1:35,'DisplayStyle','tile','ShowEmptyBins','on')
colormap(flipud(bone))
axis square
xlabel('ASA')
ylabel('neighbors')





set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



[~,all_residue_type]=structure_types(secondary_mat);


%plot ASA and neighbors for each structure and residue
v_mean_asa_structure=nan(length(structure_labels),1);
v_mean_neighbor_structure=nan(length(structure_labels),1);

for i=1:length(structure_labels)
    
    temp_idx=all_residue_type==i;
    
    v_mean_asa_structure(i)=mean(asa_mat(temp_idx),'omitnan');
    v_mean_neighbor_structure(i)=mean(neighbor_mat(temp_idx),'omitnan');
    
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
bar(v_mean_asa_structure)
title('ASA')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
axis square



subplot(2,4,2)
bar(v_mean_neighbor_structure)
title('neighbors')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
axis square



v_mean_asa_residue=nan(length(aa_labels),1);
v_mean_neighbor_residue=nan(length(aa_labels),1);
for i=1:length(aa_labels)
    
    temp_idx=ismember(residue_mat,aa_labels{i});
    
    v_mean_asa_residue(i)=mean(asa_mat(temp_idx),'omitnan');
    v_mean_neighbor_residue(i)=mean(neighbor_mat(temp_idx),'omitnan');
    
end

subplot(2,2,3)
bar(v_mean_asa_residue)
title('ASA')
xticks(1:length(aa_labels))
xticklabels(aa_labels)

subplot(2,2,4)
bar(v_mean_neighbor_residue)
title('neighbors')
xticks(1:length(aa_labels))
xticklabels(aa_labels)





set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;




figure('units','normalized','outerposition',[0 0 1 1])

%break down the scatter by secondary structure
for i=1:length(structure_labels)
    
    temp_idx=all_residue_type==i;
    
    v1=reshape(asa_mat(temp_idx),1,[]);
    v1=v1(~isnan(v1));
    v2=reshape(neighbor_mat(temp_idx),1,[]);
    v2=v2(~isnan(v2));
    
    subplot(2,4,i)
    %scatter(v1,v2,10,'k','filled')
    histogram2(v1,v2,0:5:250,4:1:35,'DisplayStyle','tile','ShowEmptyBins','on')
    colormap(flipud(bone))
    axis square
    xlabel('ASA')
    ylabel('neighbors')
    title(structure_labels{i})
    ylim([4 35])
    xlim([0 250])
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





%also by residue


figure('units','normalized','outerposition',[0 0 1 1])


for i=1:length(aa_labels)
    
    temp_idx=ismember(residue_mat,aa_labels{i});
    
    v1=reshape(asa_mat(temp_idx),1,[]);
    v1=v1(~isnan(v1));
    v2=reshape(neighbor_mat(temp_idx),1,[]);
    v2=v2(~isnan(v2));
    
    subplot(3,7,i)
    %scatter(v1,v2,10,'k','filled')
    histogram2(v1,v2,0:5:250,4:1:35,'DisplayStyle','tile','ShowEmptyBins','on')
    colormap(flipud(bone))
    axis square
    xlabel('ASA')
    ylabel('neighbors')
    title(aa_labels{i})
    ylim([4 35])
    xlim([0 250])
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;








%collect global stats and break down by secondary/residue
gene_structure_asa_mat_1K=nan(length(genes_to_use),length(structure_labels));
gene_structure_asa_mat_sim=nan(length(genes_to_use),length(structure_labels));

gene_structure_neighbor_mat_1K=nan(length(genes_to_use),length(structure_labels));
gene_structure_neighbor_mat_sim=nan(length(genes_to_use),length(structure_labels));


gene_residue_asa_mat_1K=nan(length(genes_to_use),length(aa_labels));
gene_residue_asa_mat_sim=nan(length(genes_to_use),length(aa_labels));

gene_residue_neighbor_mat_1K=nan(length(genes_to_use),length(aa_labels));
gene_residue_neighbor_mat_sim=nan(length(genes_to_use),length(aa_labels));

for i=1:length(genes_to_use)

    cols_to_use=~isnan(asa_mat(i,:));
    protein_length(i)=sum(cols_to_use);

    temp_asa=asa_mat(i,cols_to_use);
    temp_neighbor=neighbor_mat(i,cols_to_use);
    temp_residue=residue_mat(i,cols_to_use);
    [~,temp_secondary]=structure_types(secondary_mat(i,cols_to_use));
    
    pos_1K=mutation_mat_1K(i,:);
    pos_1K=pos_1K(~isnan(pos_1K));
    pos_1K(pos_1K>protein_length(i))=[];

    pos_sim=mutation_mat_simulated(i,:);
    pos_sim=pos_sim(~isnan(pos_sim));
    
    
    temp_asa_1K=temp_asa(pos_1K);
    temp_asa_sim=temp_asa(pos_sim);
    
    temp_neighbor_1K=temp_neighbor(pos_1K);
    temp_neighbor_sim=temp_neighbor(pos_sim);
    
    temp_secondary_1K=temp_secondary(pos_1K);
    temp_secondary_sim=temp_secondary(pos_sim);
    
    temp_residue_1K=temp_residue(pos_1K);
    temp_residue_sim=temp_residue(pos_sim);
    %mean asa by gene and by residue/secondary structure
    for j=1:length(structure_labels)
        
        temp_structure_idx=temp_secondary_1K==j;
        gene_structure_asa_mat_1K(i,j)=mean(temp_asa_1K(temp_structure_idx));
        gene_structure_neighbor_mat_1K(i,j)=mean(temp_neighbor_1K(temp_structure_idx));
        
        temp_structure_idx=temp_secondary_sim==j;
        gene_structure_asa_mat_sim(i,j)=mean(temp_asa_sim(temp_structure_idx));
        gene_structure_neighbor_mat_sim(i,j)=mean(temp_neighbor_sim(temp_structure_idx));
        
    end
    
    for j=1:length(aa_labels)
        
        temp_residue_idx=ismember(temp_residue_1K,aa_labels{j});
        gene_residue_asa_mat_1K(i,j)=mean(temp_asa_1K(temp_residue_idx));
        gene_residue_neighbor_mat_1K(i,j)=mean(temp_neighbor_1K(temp_residue_idx));
        
        temp_residue_idx=ismember(temp_residue_sim,aa_labels{j});
        gene_residue_asa_mat_sim(i,j)=mean(temp_asa_sim(temp_residue_idx));
        gene_residue_neighbor_mat_sim(i,j)=mean(temp_neighbor_sim(temp_residue_idx));
        
    end

end




mean_structure_asa_relative=mean(gene_structure_asa_mat_1K./...
    gene_structure_asa_mat_sim,'omitnan');

mean_structure_neighbor_relative=mean(gene_structure_neighbor_mat_1K./...
    gene_structure_neighbor_mat_sim,'omitnan');



mean_residue_asa_relative=mean(gene_residue_asa_mat_1K./...
    gene_residue_asa_mat_sim,'omitnan');

mean_residue_neighbor_relative=mean(gene_residue_neighbor_mat_1K./...
    gene_residue_neighbor_mat_sim,'omitnan');




%plot basic summary statistics for all residues/secondary structures first


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
bar(mean_structure_asa_relative,'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('ASA')
[h p]=ttest(gene_structure_asa_mat_1K./...
    gene_structure_asa_mat_sim-1);
for i=1:length(structure_labels)
    text(i,0.85,num2str(p(i)))
end

subplot(2,2,2)
hold on
bar(mean_structure_neighbor_relative,'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('neighbors')
[h p]=ttest(gene_structure_neighbor_mat_1K./...
    gene_structure_neighbor_mat_sim-1);
for i=1:length(structure_labels)
    text(i,0.85,num2str(p(i)))
end
 


subplot(2,2,3)
hold on
bar(mean_residue_asa_relative,'BaseValue',1)
ylim([0.5 1.5])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('ASA')
[h p]=ttest(gene_residue_asa_mat_1K./...
    gene_residue_asa_mat_sim-1);
for i=1:length(aa_labels)
    text(i,0.85,num2str(p(i)))
end



subplot(2,2,4)
hold on
bar(mean_residue_neighbor_relative,'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('neighbors')
[h p]=ttest(gene_residue_neighbor_mat_1K./...
    gene_residue_neighbor_mat_sim-1);
for i=1:length(aa_labels)
    text(i,0.85,num2str(p(i)))
end




set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





%break down by essential/haploinsufficient/etc and by gene age
ess_idx=ismember(genes_to_use,essential_table.GeneSystematicName);

hap_idx=ismember(genes_to_use,hap_table.ORF);

%ess_idx=hap_idx;


mean_structure_asa_relative_ess=mean(gene_structure_asa_mat_1K(ess_idx,:)./...
    gene_structure_asa_mat_sim(ess_idx,:),'omitnan');
mean_structure_asa_relative_non_ess=mean(gene_structure_asa_mat_1K(~ess_idx,:)./...
    gene_structure_asa_mat_sim(~ess_idx,:),'omitnan');

mean_structure_neighbor_relative_ess=mean(gene_structure_neighbor_mat_1K(ess_idx,:)./...
    gene_structure_neighbor_mat_sim(ess_idx,:),'omitnan');
mean_structure_neighbor_relative_non_ess=mean(gene_structure_neighbor_mat_1K(~ess_idx,:)./...
    gene_structure_neighbor_mat_sim(~ess_idx,:),'omitnan');



mean_residue_asa_relative_ess=mean(gene_residue_asa_mat_1K(ess_idx,:)./...
    gene_residue_asa_mat_sim(ess_idx,:),'omitnan');
mean_residue_asa_relative_non_ess=mean(gene_residue_asa_mat_1K(~ess_idx,:)./...
    gene_residue_asa_mat_sim(~ess_idx,:),'omitnan');

mean_residue_neighbor_relative_ess=mean(gene_residue_neighbor_mat_1K(ess_idx,:)./...
    gene_residue_neighbor_mat_sim(ess_idx,:),'omitnan');
mean_residue_neighbor_relative_non_ess=mean(gene_residue_neighbor_mat_1K(~ess_idx,:)./...
    gene_residue_neighbor_mat_sim(~ess_idx,:),'omitnan');


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
bar([mean_structure_asa_relative_ess' mean_structure_asa_relative_non_ess'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('ASA')
legend({'essential','other'})




subplot(2,2,2)
hold on
bar([mean_structure_neighbor_relative_ess' mean_structure_neighbor_relative_non_ess'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('neighbors')
legend({'essential','other'})
 




subplot(2,2,3)
hold on
bar([mean_residue_asa_relative_ess' mean_residue_asa_relative_non_ess'],'BaseValue',1)
ylim([0.5 1.5])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('ASA')
legend({'essential','other'})



subplot(2,2,4)
hold on
bar([mean_residue_neighbor_relative_ess' mean_residue_neighbor_relative_non_ess'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('neighbors')
legend({'essential','other'})




set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



gene_age_info=readtable([dependency_directory 'Supplementary_Data_4_Doughty_et_al_2020.xlsx']);

group_labels={'Group I','Group II','Group III','WGD','Group IV','Group V'};

gene_age=nan(length(genes_to_use),1);
for i=1:length(genes_to_use)

    temp_idx=ismember(gene_age_info.Var1,genes_to_use{i});

    if sum(temp_idx)>0

        gene_age(i)=find(ismember(group_labels,gene_age_info.Var2{temp_idx}));

    end

end

young_idx=gene_age==6;
old_idx=gene_age==1;

mean_structure_asa_relative_old=mean(gene_structure_asa_mat_1K(old_idx,:)./...
    gene_structure_asa_mat_sim(old_idx,:),'omitnan');
mean_structure_asa_relative_young=mean(gene_structure_asa_mat_1K(young_idx,:)./...
    gene_structure_asa_mat_sim(young_idx,:),'omitnan');

mean_structure_neighbor_relative_old=mean(gene_structure_neighbor_mat_1K(old_idx,:)./...
    gene_structure_neighbor_mat_sim(old_idx,:),'omitnan');
mean_structure_neighbor_relative_young=mean(gene_structure_neighbor_mat_1K(young_idx,:)./...
    gene_structure_neighbor_mat_sim(young_idx,:),'omitnan');



mean_residue_asa_relative_old=mean(gene_residue_asa_mat_1K(old_idx,:)./...
    gene_residue_asa_mat_sim(old_idx,:),'omitnan');
mean_residue_asa_relative_young=mean(gene_residue_asa_mat_1K(young_idx,:)./...
    gene_residue_asa_mat_sim(young_idx,:),'omitnan');

mean_residue_neighbor_relative_old=mean(gene_residue_neighbor_mat_1K(old_idx,:)./...
    gene_residue_neighbor_mat_sim(old_idx,:),'omitnan');
mean_residue_neighbor_relative_young=mean(gene_residue_neighbor_mat_1K(young_idx,:)./...
    gene_residue_neighbor_mat_sim(young_idx,:),'omitnan');


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
bar([mean_structure_asa_relative_old' mean_structure_asa_relative_young'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('ASA')
legend({'old','young'})




subplot(2,2,2)
hold on
bar([mean_structure_neighbor_relative_old' mean_structure_neighbor_relative_young'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('1K/sim')
title('neighbors')
legend({'old','young'})
 




subplot(2,2,3)
hold on
bar([mean_residue_asa_relative_old' mean_residue_asa_relative_young'],'BaseValue',1)
ylim([0.5 1.5])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('ASA')
legend({'old','young'})



subplot(2,2,4)
hold on
bar([mean_residue_neighbor_relative_old' mean_residue_neighbor_relative_young'],'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('1K/sim')
title('neighbors')
legend({'old','young'})



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





%also example 1K vs sim ASA and neighbor histograms for interesting residues/structures
figure('units','normalized','outerposition',[0 0 1 1])
m=1;

%alphahelix vs betaladder vs unstr
structure_to_plot=[1 5 8];

for i=1:length(structure_to_plot)
    
    v1=gene_structure_asa_mat_1K(:,structure_to_plot(i));
    v1=v1(~isnan(v1));
    v2=gene_structure_asa_mat_sim(:,structure_to_plot(i));
    v2=v2(~isnan(v2));
    
    subplot(2,6,m)
    hold on
    histogram(v1,0:10:250,'Normalization','probability')
    histogram(v2,0:10:250,'Normalization','probability')
    axis square
    legend({'1K','sim'})
    title(structure_labels{structure_to_plot(i)})
    ylabel('relative frequency')
    xlabel('ASA')
    m=m+1;

    
    v1=gene_structure_neighbor_mat_1K(:,structure_to_plot(i));
    v1=v1(~isnan(v1));
    v2=gene_structure_neighbor_mat_sim(:,structure_to_plot(i));
    v2=v2(~isnan(v2));
    
    subplot(2,6,m)
    hold on
    histogram(v1,0:2:35,'Normalization','probability')
    histogram(v2,0:2:35,'Normalization','probability')
    axis square
    legend({'1K','sim'})
    title(structure_labels{structure_to_plot(i)})
    ylabel('relative frequency')
    xlabel('neighbors')
    m=m+1;
    
end





m=7;
aa_to_plot={'A','C','K'};

for i=1:length(aa_to_plot)
    
    aa_idx=find(ismember(aa_labels,aa_to_plot{i}));
    
    v1=gene_residue_asa_mat_1K(:,aa_idx);
    v1=v1(~isnan(v1));
    v2=gene_residue_asa_mat_sim(:,aa_idx);
    v2=v2(~isnan(v2));
    
    subplot(2,6,m)
    hold on
    histogram(v1,0:10:250,'Normalization','probability')
    histogram(v2,0:10:250,'Normalization','probability')
    axis square
    legend({'1K','sim'})
    title(aa_to_plot{i})
    ylabel('relative frequency')
    xlabel('ASA')
    m=m+1;

    
    
    v1=gene_residue_neighbor_mat_1K(:,aa_idx);
    v1=v1(~isnan(v1));
    v2=gene_residue_neighbor_mat_sim(:,aa_idx);
    v2=v2(~isnan(v2));
    
    subplot(2,6,m)
    hold on
    histogram(v1,0:2:35,'Normalization','probability')
    histogram(v2,0:2:35,'Normalization','probability')
    axis square
    legend({'1K','sim'})
    title(aa_to_plot{i})
    ylabel('relative frequency')
    xlabel('neighbors')
    m=m+1;
    
end





set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;










asa_mat_1K=nan(size(af_mat_1K));
neighbor_mat_1K=nan(size(af_mat_1K));
residue_mat_1K=cell(size(af_mat_1K));
structure_mat_1K=nan(size(af_mat_1K));

asa_mat_sim=nan(size(mutation_mat_simulated));
neighbor_mat_sim=nan(size(mutation_mat_simulated));
residue_mat_sim=cell(size(mutation_mat_simulated));
structure_mat_sim=nan(size(mutation_mat_simulated));
%gene-wise stats
for i=1:length(genes_to_use)
    
    cols_to_use=~isnan(asa_mat(i,:));
    
    temp_asa=asa_mat(i,cols_to_use);
    temp_neighbor=neighbor_mat(i,cols_to_use);
    temp_residue=residue_mat(i,cols_to_use);
    [~,temp_secondary]=structure_types(secondary_mat(i,cols_to_use));
    
    pos_1K=mutation_mat_1K(i,:);
    pos_1K=pos_1K(~isnan(pos_1K));
    pos_1K(pos_1K>protein_length(i))=[];

    pos_sim=mutation_mat_simulated(i,:);
    pos_sim=pos_sim(~isnan(pos_sim));
    
    %also build matrices to compare to AF
    asa_mat_1K(i,1:length(pos_1K))=temp_asa(pos_1K);
    asa_mat_sim(i,1:length(pos_sim))=temp_asa(pos_sim);
    
    asa_1K(i)=mean(temp_asa(pos_1K));
    asa_sim(i)=mean(temp_asa(pos_sim));
    rel_asa(i)=asa_1K(i)/asa_sim(i);
    
    
    neighbor_mat_1K(i,1:length(pos_1K))=temp_neighbor(pos_1K);
    neighbor_mat_sim(i,1:length(pos_sim))=temp_neighbor(pos_sim);
    
    neighbor_1K(i)=mean(temp_neighbor(pos_1K));
    neighbor_sim(i)=mean(temp_neighbor(pos_sim));
    rel_neighbor(i)=neighbor_1K(i)/neighbor_sim(i);
    
    residue_mat_1K(i,1:length(pos_1K))=temp_residue(pos_1K);
    structure_mat_1K(i,1:length(pos_1K))=temp_secondary(pos_1K);
        
    residue_mat_sim(i,1:length(pos_sim))=temp_residue(pos_sim);
    structure_mat_sim(i,1:length(pos_sim))=temp_secondary(pos_sim);
    
    %also how many mutations fixed relative to all possible
    %to control for general essentiality
    mutations_1K(i)=length(pos_1K)/protein_length(i);
    mutations_sim(i)=length(pos_sim)/protein_length(i);
    rel_mutations(i)=length(pos_1K)/length(pos_sim);
    
end



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
hold on
v1=asa_1K;
v1=v1(~isnan(v1));
v2=asa_sim;
v2=v2(~isnan(v2));
histogram(v1,0:10:250)%,'Normalization','probability')
histogram(v2,0:10:250)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('ASA')
ylabel('number of genes')



subplot(2,4,2)
hold on
v1=neighbor_1K;
v1=v1(~isnan(v1));
v2=neighbor_sim;
v2=v2(~isnan(v2));
histogram(v1,0:2:35)%,'Normalization','probability')
histogram(v2,0:2:35)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('neighbors')
ylabel('number of genes')


subplot(2,4,3)
hold on
v1=mutations_1K;
v1=v1(~isnan(v1));
v2=mutations_sim;
v2=v2(~isnan(v2));
histogram(v1)%,'Normalization','probability')
histogram(v2)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('mutations/residue')
ylabel('number of genes')





subplot(2,4,5)
hold on
%scatter(rel_asa,rel_neighbor,10,'k','filled')
histogram2(rel_asa,rel_neighbor,0.5:0.02:2,0.5:0.02:1.5,'DisplayStyle','tile','ShowEmptyBins','on')
colormap(flipud(bone))
xlim([0.5 2])
ylim([0.5 1.5])
axis square
xlabel('relative ASA 1K/sim')
ylabel('relative neighbors 1K/sim')


subplot(2,4,6)
hold on
%scatter(rel_mutations,rel_asa,10,'k','filled')
histogram2(rel_mutations,rel_asa,0:0.002:0.1,0.5:0.02:2,'DisplayStyle','tile','ShowEmptyBins','on')
colormap(flipud(bone))
ylim([0.5 2])
xlim([0 0.1])
axis square
xlabel('relative # mutations 1K/sim')
ylabel('relative ASA 1K/sim')


subplot(2,4,7)
hold on
%scatter(rel_mutations,rel_neighbor,10,'k','filled')
histogram2(rel_mutations,rel_neighbor,0:0.002:0.1,0.5:0.02:1.5,'DisplayStyle','tile','ShowEmptyBins','on')
colormap(flipud(bone))
ylim([0.5 1.5])
xlim([0 0.1])
axis square
xlabel('relative # mutations 1K/sim')
ylabel('relative neighbors 1K/sim')





set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;






n_to_output=500;
%identify most and least constrained genes
[v_sorted,sort_idx]=sort(rel_mutations,'ascend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'most_constrained_genes_mutations.txt'])


[v_sorted,sort_idx]=sort(rel_mutations,'descend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'least_constrained_genes_mutations.txt'])




[v_sorted,sort_idx]=sort(rel_asa,'ascend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'least_constrained_genes_asa.txt'])


[v_sorted,sort_idx]=sort(rel_asa,'descend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'most_constrained_genes_asa.txt'])




[v_sorted,sort_idx]=sort(rel_neighbor,'descend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'least_constrained_genes_neighbor.txt'])


[v_sorted,sort_idx]=sort(rel_neighbor,'ascend');

to_output=table(genes_to_use(sort_idx(1:n_to_output)));
writetable(to_output,[output_directory 'most_constrained_genes_neighbor.txt'])


%close all

residue_mat_1K(cellfun(@isempty,residue_mat_1K))={'NA'};
residue_mat_sim(cellfun(@isempty,residue_mat_sim))={'NA'};

[~,structure_mat]=structure_types(secondary_mat);

%2D analysis by residue/structure
mean_2D_asa_relative=nan(length(aa_labels),length(structure_labels));
mean_2D_neighbor_relative=nan(length(aa_labels),length(structure_labels));
for i=1:length(aa_labels)
    
    temp_residue_idx1=ismember(residue_mat_1K,aa_labels{i});
    temp_residue_idx2=ismember(residue_mat_sim,aa_labels{i});
    
    for j=1:length(structure_labels)
        
        temp_structure_idx1=structure_mat_1K==j;
        
        temp_idx_to_use1=logical(temp_residue_idx1.*temp_structure_idx1);
        
        
        temp_structure_idx2=structure_mat_sim==j;
        
        temp_idx_to_use2=logical(temp_residue_idx2.*temp_structure_idx2);
        
        
        mean_2D_asa_relative(i,j)=mean(asa_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(asa_mat_sim(temp_idx_to_use2),'omitnan');
        mean_2D_neighbor_relative(i,j)=mean(neighbor_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(neighbor_mat_sim(temp_idx_to_use2),'omitnan');
        
    end
    
end



figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,6,1)
imagesc(mean_2D_asa_relative,[1 1.7])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('ASA 1K/sim')

m = size(get(gcf,'colormap'),1);
%red to blue colormap
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
%colormap(flipud(c))
colormap(c)
colorbar


subplot(2,6,2)
imagesc(mean_2D_neighbor_relative,[0.7 1])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
colormap(c)
colorbar
title('neighbors 1K/sim')


subplot(2,4,3)
v1=reshape(mean_2D_neighbor_relative,[],1);
v2=reshape(mean_2D_asa_relative,[],1);
scatter(v1,v2,25,'k','filled');
axis square
xlabel('neighbors 1K/sim')
ylabel('ASA 1K/sim')
xlim([0.7 1])
ylim([1 1.7])






set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





[n_proteins,max_residues]=size(asa_mat);
[~,n_cols]=size(domain_mat);
domain_mat(:,(n_cols+1):max_residues)=0;

[n_rows,~]=size(domain_mat);
domain_mat((n_rows+1):n_proteins,:)=0;

%also integrate information on whether residue is within SuperFam domain
%gene-wise stats
domain_mat_1K=nan(size(asa_mat_1K));
domain_mat_sim=nan(size(asa_mat_sim));
for i=1:length(genes_to_use)
    
    cols_to_use=~isnan(asa_mat(i,:));
    
    temp_asa=asa_mat(i,cols_to_use);
    temp_neighbor=neighbor_mat(i,cols_to_use);
    temp_residue=residue_mat(i,cols_to_use);
    [~,temp_secondary]=structure_types(secondary_mat(i,cols_to_use));
    
    v_domain=logical(domain_mat(i,find(cols_to_use)));
    
    pos_1K=mutation_mat_1K(i,:);
    pos_1K=pos_1K(~isnan(pos_1K));
    pos_1K(pos_1K>protein_length(i))=[];
    
    v_domain_1K=logical(domain_mat(i,pos_1K));
    domain_mat_1K(i,1:length(v_domain_1K))=v_domain_1K;

    pos_sim=mutation_mat_simulated(i,:);
    pos_sim=pos_sim(~isnan(pos_sim));
    
    v_domain_sim=logical(domain_mat(i,pos_sim));
    domain_mat_sim(i,1:length(v_domain_sim))=v_domain_sim;
    
    temp_asa_1K=temp_asa(pos_1K);
    
    asa_1K_domain(i)=mean(temp_asa_1K(v_domain_1K));
    asa_1K_nondomain(i)=mean(temp_asa_1K(~v_domain_1K));
    
    temp_asa_sim=temp_asa(pos_sim);
    
    asa_sim_domain(i)=mean(temp_asa_sim(v_domain_sim));
    asa_sim_nondomain(i)=mean(temp_asa_sim(~v_domain_sim));
    
    
    rel_asa_domain(i)=asa_1K_domain(i)/asa_sim_domain(i);
    rel_asa_nondomain(i)=asa_1K_nondomain(i)/asa_sim_nondomain(i);
    
    
    
    temp_neighbor_1K=temp_neighbor(pos_1K);
    
    neighbor_1K_domain(i)=mean(temp_neighbor_1K(v_domain_1K));
    neighbor_1K_nondomain(i)=mean(temp_neighbor_1K(~v_domain_1K));
    
    temp_neighbor_sim=temp_neighbor(pos_sim);
    
    neighbor_sim_domain(i)=mean(temp_neighbor_sim(v_domain_sim));
    neighbor_sim_nondomain(i)=mean(temp_neighbor_sim(~v_domain_sim));
    
    
    rel_neighbor_domain(i)=neighbor_1K_domain(i)/neighbor_sim_domain(i);
    rel_neighbor_nondomain(i)=neighbor_1K_nondomain(i)/neighbor_sim_nondomain(i);
    
    
    %also how many mutations fixed relative to all possible
    %to control for general essentiality
    
    
    
    mutations_1K_domain(i)=sum(v_domain_1K)/sum(v_domain);
    mutations_1K_nondomain(i)=sum(~v_domain_1K)/sum(~v_domain);
    
    mutations_sim_domain(i)=sum(v_domain_sim)/sum(v_domain);
    mutations_sim_nondomain(i)=sum(~v_domain_sim)/sum(~v_domain);
    
    rel_mutations_domain(i)=mutations_1K_domain(i)/mutations_sim_domain(i);
    rel_mutations_nondomain(i)=mutations_1K_nondomain(i)/mutations_sim_nondomain(i);
    
end




figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,4,1)
hold on
v1=rel_asa_domain;
v1(isnan(v1))=[];
to_plot{1}=v1;

v2=rel_asa_nondomain;
v2(isnan(v2))=[];
to_plot{2}=v2;

histogram(v1,0:0.05:2,'Normalization','probability')
histogram(v2,0:0.05:2,'Normalization','probability')
legend({'in domain','outside domain'})
title('ASA')
axis square
xlabel('1K/sim')


subplot(2,4,2)
hold on
v1=rel_neighbor_domain;
v1(isnan(v1))=[];
to_plot{3}=v1;

v2=rel_neighbor_nondomain;
v2(isnan(v2))=[];
to_plot{4}=v2;

histogram(v1,0:0.05:2,'Normalization','probability')
histogram(v2,0:0.05:2,'Normalization','probability')
legend({'in domain','outside domain'})
title('neighbors')
axis square
xlabel('1K/sim')



subplot(2,4,3)
hold on
easy_box(to_plot)
axis square
ylim([0 2])
plot(xlim,[1 1],':r')
for i=1:length(to_plot)
    [h p]=ttest(to_plot{i}-1);
    text(i,0.2,num2str(p))
end
xtickangle(45)
xticklabels({'ASA in domain','ASA outside',...
    'neighbor in domain','neighbor outside'})
xlabel('1K/sim')



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



%rerun per residue and per structure analysis for in domain vs not

%2D analysis by residue/structure
mean_2D_asa_relative_domain=nan(length(aa_labels),length(structure_labels));
mean_2D_asa_relative_nondomain=nan(length(aa_labels),length(structure_labels));

mean_2D_neighbor_relative_domain=nan(length(aa_labels),length(structure_labels));
mean_2D_neighbor_relative_nondomain=nan(length(aa_labels),length(structure_labels));

domain_mat_1K(isnan(domain_mat_1K))=0;
domain_mat_sim(isnan(domain_mat_sim))=0;
for i=1:length(aa_labels)
    
    for j=1:length(structure_labels)
        
        temp_residue_idx=ismember(residue_mat_1K,aa_labels{i});
        temp_structure_idx=structure_mat_1K==j;
        
        temp_idx_to_use1=logical(temp_residue_idx.*temp_structure_idx.*domain_mat_1K);
        
        
        temp_residue_idx=ismember(residue_mat_sim,aa_labels{i});
        temp_structure_idx=structure_mat_sim==j;
        
        temp_idx_to_use2=logical(temp_residue_idx.*temp_structure_idx.*domain_mat_sim);
        
        
        mean_2D_asa_relative_domain(i,j)=mean(asa_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(asa_mat_sim(temp_idx_to_use2),'omitnan');
        mean_2D_neighbor_relative_domain(i,j)=mean(neighbor_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(neighbor_mat_sim(temp_idx_to_use2),'omitnan');
        
        
        
        temp_residue_idx=ismember(residue_mat_1K,aa_labels{i});
        temp_structure_idx=structure_mat_1K==j;
        
        temp_idx_to_use1=logical(temp_residue_idx.*temp_structure_idx.*~domain_mat_1K);
        
        
        temp_residue_idx=ismember(residue_mat_sim,aa_labels{i});
        temp_structure_idx=structure_mat_sim==j;
        
        temp_idx_to_use2=logical(temp_residue_idx.*temp_structure_idx.*~domain_mat_sim);
        
        
        mean_2D_asa_relative_nondomain(i,j)=mean(asa_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(asa_mat_sim(temp_idx_to_use2),'omitnan');
        mean_2D_neighbor_relative_nondomain(i,j)=mean(neighbor_mat_1K(temp_idx_to_use1),'omitnan')/...
            mean(neighbor_mat_sim(temp_idx_to_use2),'omitnan');
        
    end
    
end






figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,6,1)
imagesc(mean_2D_asa_relative_domain,[1 1.7])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('ASA 1K/sim in domain')
colormap(c)
colorbar

subplot(2,6,2)
imagesc(mean_2D_asa_relative_nondomain,[1 1.7])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('ASA 1K/sim outside domain')
colormap(c)
colorbar




subplot(2,6,3)
imagesc(mean_2D_neighbor_relative_domain,[0.7 1])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
colormap(c)
colorbar
title('neighbors 1K/sim in domain')



subplot(2,6,4)
imagesc(mean_2D_neighbor_relative_nondomain,[0.7 1])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
colormap(c)
colorbar
title('neighbors 1K/sim outside domain')





subplot(2,4,5)
v1=reshape(mean_2D_neighbor_relative_domain,[],1);
v2=reshape(mean_2D_asa_relative_domain,[],1);
scatter(v1,v2,25,'k','filled');
axis square
xlabel('neighbors 1K/sim')
ylabel('ASA 1K/sim')
title('in domain')
xlim([0.7 1])
ylim([1 1.7])


subplot(2,4,6)
v1=reshape(mean_2D_neighbor_relative_nondomain,[],1);
v2=reshape(mean_2D_asa_relative_nondomain,[],1);
scatter(v1,v2,25,'k','filled');
axis square
xlabel('neighbors 1K/sim')
ylabel('ASA 1K/sim')
title('outside domain')
xlim([0.7 1])
ylim([1 1.7])






set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;








%start AF analysis
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
v1=af_mat_1K;
v1=v1(~isnan(v1));
histogram(v1)
%set(gca,'XScale','log')
set(gca,'YScale','log')
axis square
xlabel('allele frequency in 1K')
ylabel('number of missense variants')


%number of variants per gene
v1=sum(~isnan(af_mat_1K),2);
subplot(2,4,2)
histogram(v1)
axis square
xlabel('missense variants per gene')
ylabel('number of genes')


subplot(2,4,3)
hold on
%bin by AF
af_bins=[0 1e-3 1e-2 Inf];

n_variants=nan(length(af_bins)-1,1);
clear to_plot
for i=1:(length(af_bins)-1)
    
    temp_idx=logical((af_mat_1K>af_bins(i)).*(af_mat_1K<=af_bins(i+1)));
    n_variants(i)=sum(sum(temp_idx));
    
    to_plot{i}=asa_mat_1K(temp_idx);
    to_plot{i}=to_plot{i}(~isnan(to_plot{i}));
    
    histogram(to_plot{i},0:10:250,'Normalization','probability')
    
end
%set(gca,'YScale','log')
axis square
xlabel('ASA')
ylabel('rel. freq.')
legend({'AF<0.1%','0.1%<AF<1%','AF>1%'})




subplot(2,4,4)
hold on
n_variants=nan(length(af_bins)-1,1);
clear to_plot
for i=1:(length(af_bins)-1)
    
    temp_idx=logical((af_mat_1K>af_bins(i)).*(af_mat_1K<=af_bins(i+1)));
    n_variants(i)=sum(sum(temp_idx));
    
    to_plot{i}=neighbor_mat_1K(temp_idx);
    to_plot{i}=to_plot{i}(~isnan(to_plot{i}));
    
    histogram(to_plot{i},0:2:35,'Normalization','probability')
    
end
%set(gca,'YScale','log')
axis square
xlabel('neighbors')
ylabel('rel. freq.')
legend({'AF<0.1%','0.1%<AF<1%','AF>1%'})




%calibrate AF threshold
%p value and effect size
thresh_to_test=10.^(-5:0.1:-1);

asa_p_val=nan(length(thresh_to_test),1);
asa_effect_size=nan(length(thresh_to_test),1);
for i=1:length(thresh_to_test)
    
    temp_idx=af_mat_1K<=thresh_to_test(i);
    
    if sum(sum(temp_idx))>0
    
        [h p]=ttest2(asa_mat_1K(temp_idx),asa_mat_1K(~temp_idx));

        asa_p_val(i)=p;
        asa_effect_size(i)=mean(asa_mat_1K(temp_idx),'omitnan')./...
            mean(asa_mat_1K(~temp_idx),'omitnan');
        
        
    
        [h p]=ttest2(neighbor_mat_1K(temp_idx),neighbor_mat_1K(~temp_idx));

        neighbor_p_val(i)=p;
        neighbor_effect_size(i)=mean(neighbor_mat_1K(temp_idx),'omitnan')./...
            mean(neighbor_mat_1K(~temp_idx),'omitnan');
        
    end
    
end

subplot(2,4,5)
hold on
scatter(asa_effect_size,-log10(asa_p_val),10,'k','filled')
xlim([0.8 1])
xlabel('effect size (rare/common)')
ylabel('-log10 p')
title('ASA')
ylim([0 350])
axis square


subplot(2,4,6)
hold on
scatter(neighbor_effect_size,-log10(neighbor_p_val),10,'k','filled')
xlim([1 1.2])
xlabel('effect size (rare/common)')
ylabel('-log10 p')
title('neighbors')
ylim([0 350])
axis square




set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





figure('units','normalized','outerposition',[0 0 1 1])

asa_p_val=nan(length(structure_labels),length(thresh_to_test));
asa_effect_size=nan(length(structure_labels),length(thresh_to_test));
neighbor_p_val=nan(length(structure_labels),length(thresh_to_test));
neighbor_effect_size=nan(length(structure_labels),length(thresh_to_test));
%repeat this analysis subsetting by structure
for i=1:length(structure_labels)
    
    temp_idx1=structure_mat_1K==i;
    
    for j=1:length(thresh_to_test)
    
        temp_idx2=af_mat_1K<=thresh_to_test(j);
        
        temp_idx3=logical(temp_idx1.*temp_idx2);
        temp_idx4=logical(temp_idx1.*~temp_idx2);
        
        if sum(sum(temp_idx))>0
    
            [h p]=ttest2(asa_mat_1K(temp_idx3),asa_mat_1K(temp_idx4));

            asa_p_val(i,j)=p;
            asa_effect_size(i,j)=mean(asa_mat_1K(temp_idx3),'omitnan')./...
                mean(asa_mat_1K(temp_idx4),'omitnan');



            [h p]=ttest2(neighbor_mat_1K(temp_idx3),neighbor_mat_1K(temp_idx4));

            neighbor_p_val(i,j)=p;
            neighbor_effect_size(i,j)=mean(neighbor_mat_1K(temp_idx3),'omitnan')./...
                mean(neighbor_mat_1K(temp_idx4),'omitnan');

        end
        
    end
    
    subplot(2,4,i)
    hold on
    scatter(asa_effect_size(i,:),-log10(asa_p_val(i,:)),50,'.k')
    scatter(neighbor_effect_size(i,:),-log10(neighbor_p_val(i,:)),25,'ok')
    xlim([0.85 1.15])
    xlabel('effect size (rare/common)')
    ylabel('-log10 p')
    ylim([0 Inf])
    title(structure_labels{i})
    legend({'ASA','neighbors'})
    axis square
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



%also for residues



figure('units','normalized','outerposition',[0 0 1 1])

asa_p_val=nan(length(aa_labels),length(thresh_to_test));
asa_effect_size=nan(length(aa_labels),length(thresh_to_test));
neighbor_p_val=nan(length(aa_labels),length(thresh_to_test));
neighbor_effect_size=nan(length(aa_labels),length(thresh_to_test));

for i=1:length(aa_labels)
    
    temp_idx1=ismember(residue_mat_1K,aa_labels{i});
    
    for j=1:length(thresh_to_test)
    
        temp_idx2=af_mat_1K<=thresh_to_test(j);
        
        temp_idx3=logical(temp_idx1.*temp_idx2);
        temp_idx4=logical(temp_idx1.*~temp_idx2);
        
        if sum(sum(temp_idx))>0
    
            [h p]=ttest2(asa_mat_1K(temp_idx3),asa_mat_1K(temp_idx4));

            asa_p_val(i,j)=p;
            asa_effect_size(i,j)=mean(asa_mat_1K(temp_idx3),'omitnan')./...
                mean(asa_mat_1K(temp_idx4),'omitnan');



            [h p]=ttest2(neighbor_mat_1K(temp_idx3),neighbor_mat_1K(temp_idx4));

            neighbor_p_val(i,j)=p;
            neighbor_effect_size(i,j)=mean(neighbor_mat_1K(temp_idx3),'omitnan')./...
                mean(neighbor_mat_1K(temp_idx4),'omitnan');

        end
        
    end
    
    subplot(3,7,i)
    hold on
    scatter(asa_effect_size(i,:),-log10(asa_p_val(i,:)),50,'.k')
    scatter(neighbor_effect_size(i,:),-log10(neighbor_p_val(i,:)),25,'ok')
    xlim([0.7 1.3])
    xlabel('effect size (rare/common)')
    ylabel('-log10 p')
    ylim([0 Inf])
    title(aa_labels{i})
    legend({'ASA','neighbors'})
    axis square
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;


%pick a threshold and make bar plots again
af_thresh=0.001;


figure('units','normalized','outerposition',[0 0 1 1])

asa_p_val=nan(length(structure_labels),1);
asa_effect_size=nan(length(structure_labels),1);
neighbor_p_val=nan(length(structure_labels),1);
neighbor_effect_size=nan(length(structure_labels),1);
%repeat this analysis subsetting by structure
for i=1:length(structure_labels)
    
    temp_idx1=structure_mat_1K==i;
    
    temp_idx2=af_mat_1K<=af_thresh;

    temp_idx3=logical(temp_idx1.*temp_idx2);
    temp_idx4=logical(temp_idx1.*~temp_idx2);
    
    [h p]=ttest2(asa_mat_1K(temp_idx3),asa_mat_1K(temp_idx4));

    asa_p_val(i)=p;
    asa_effect_size(i)=mean(asa_mat_1K(temp_idx3),'omitnan')./...
        mean(asa_mat_1K(temp_idx4),'omitnan');



    [h p]=ttest2(neighbor_mat_1K(temp_idx3),neighbor_mat_1K(temp_idx4));

    neighbor_p_val(i)=p;
    neighbor_effect_size(i)=mean(neighbor_mat_1K(temp_idx3),'omitnan')./...
        mean(neighbor_mat_1K(temp_idx4),'omitnan');

end


subplot(2,2,1)
hold on
bar(asa_effect_size,'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('rare/common')
title('ASA')




subplot(2,2,2)
hold on
bar(neighbor_effect_size,'BaseValue',1)
ylim([0.8 1.2])
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylabel('rare/common')
title('neighbors')





asa_p_val=nan(length(aa_labels),1);
asa_effect_size=nan(length(aa_labels),1);
neighbor_p_val=nan(length(aa_labels),1);
neighbor_effect_size=nan(length(aa_labels),1);
%repeat this analysis subsetting by structure
for i=1:length(aa_labels)
    
    temp_idx1=ismember(residue_mat_1K,aa_labels{i});
    
    temp_idx2=af_mat_1K<=af_thresh;

    temp_idx3=logical(temp_idx1.*temp_idx2);
    temp_idx4=logical(temp_idx1.*~temp_idx2);
    
    [h p]=ttest2(asa_mat_1K(temp_idx3),asa_mat_1K(temp_idx4));

    asa_p_val(i)=p;
    asa_effect_size(i)=mean(asa_mat_1K(temp_idx3),'omitnan')./...
        mean(asa_mat_1K(temp_idx4),'omitnan');



    [h p]=ttest2(neighbor_mat_1K(temp_idx3),neighbor_mat_1K(temp_idx4));

    neighbor_p_val(i)=p;
    neighbor_effect_size(i)=mean(neighbor_mat_1K(temp_idx3),'omitnan')./...
        mean(neighbor_mat_1K(temp_idx4),'omitnan');

end


subplot(2,2,3)
hold on
bar(asa_effect_size,'BaseValue',1)
ylim([0.7 1.3])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('rare/common')
title('ASA')




subplot(2,2,4)
hold on
bar(neighbor_effect_size,'BaseValue',1)
ylim([0.7 1.3])
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylabel('rare/common')
title('neighbors')







set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





%do old/young split again?







%close all


at_data=readtable([dependency_directory 'arabidopsis-data/1001misAnnotated.csv']);

%convert allele counts to MAF
at_af=at_data.AC;
at_af=at_af./max(at_af);
at_af(at_af>0.5)=1-at_af(at_af>0.5);

at_asa=at_data.sasa;
at_neighbor=at_data.neighbors;

[~,at_structure]=structure_types(at_data.secondary);
at_residue=at_data.residue;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
histogram(at_af,0:0.01:0.5)
set(gca,'YScale','log')
xlabel('A.t. AF')
ylabel('frequency')
axis square



at_asa_p_val=nan(length(thresh_to_test),1);
at_asa_effect_size=nan(length(thresh_to_test),1);

at_neighbor_p_val=nan(length(thresh_to_test),1);
at_neighbor_effect_size=nan(length(thresh_to_test),1);
for i=1:length(thresh_to_test)
    
    temp_idx=at_af<=thresh_to_test(i);
    
    if sum(sum(temp_idx))>0
    
        [h p]=ttest2(at_asa(temp_idx),at_asa(~temp_idx));

        at_asa_p_val(i)=p;
        at_asa_effect_size(i)=mean(at_asa(temp_idx),'omitnan')./...
            mean(at_asa(~temp_idx),'omitnan');
        
        
    
        [h p]=ttest2(at_neighbor(temp_idx),at_neighbor(~temp_idx));

        at_neighbor_p_val(i)=p;
        at_neighbor_effect_size(i)=mean(at_neighbor(temp_idx),'omitnan')./...
            mean(at_neighbor(~temp_idx),'omitnan');
        
    end
    
end


subplot(2,4,2)
hold on
scatter(at_asa_effect_size,-log10(at_asa_p_val),50,'.k')
scatter(at_neighbor_effect_size,-log10(at_neighbor_p_val),25,'ok')
xlim([0.7 1.3])
xlabel('effect size (rare/common)')
ylabel('-log10 p')
ylim([0 350])
legend({'ASA','neighbors'})
axis square
v1=at_asa_effect_size;
v2=-log10(at_asa_p_val);
for i=1:length(thresh_to_test)
    text(v1(i),v2(i),num2str(thresh_to_test(i)))
end


%0.1% looks good

at_af_thresh=0.001;
temp_idx1=at_af<=at_af_thresh;

for i=1:length(structure_labels)
    
    temp_idx2=at_structure==i;
    
    v1=at_asa(logical(temp_idx1.*temp_idx2));
    v2=at_asa(logical(~temp_idx1.*temp_idx2));
    
    at_mean_structure_asa_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    at_p_val_structure_asa_relative(i)=p;
    
    
    v1=at_neighbor(logical(temp_idx1.*temp_idx2));
    v2=at_neighbor(logical(~temp_idx1.*temp_idx2));
    
    at_mean_structure_neighbor_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    at_p_val_structure_neighbor_relative(i)=p;

end

subplot(2,4,3)
hold on
bar(at_mean_structure_asa_relative,'BaseValue',1)
title('ASA')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.8 1.2])
for i=1:length(at_p_val_structure_asa_relative)
    text(i,0.85,num2str(at_p_val_structure_asa_relative(i)))
end
axis square


subplot(2,4,4)
hold on
bar(at_mean_structure_neighbor_relative,'BaseValue',1)
title('neighbors')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.9 1.1])
for i=1:length(at_p_val_structure_neighbor_relative)
    text(i,0.92,num2str(at_p_val_structure_neighbor_relative(i)))
end
axis square




for i=1:length(aa_labels)
    
    temp_idx2=ismember(at_residue,aa_labels{i});
    
    v1=at_asa(logical(temp_idx1.*temp_idx2));
    v2=at_asa(logical(~temp_idx1.*temp_idx2));
    
    at_mean_residue_asa_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    at_p_val_residue_asa_relative(i)=p;
    
    
    v1=at_neighbor(logical(temp_idx1.*temp_idx2));
    v2=at_neighbor(logical(~temp_idx1.*temp_idx2));
    
    at_mean_residue_neighbor_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    at_p_val_residue_neighbor_relative(i)=p;

end

subplot(2,2,3)
hold on
bar(at_mean_residue_asa_relative,'BaseValue',1)
title('ASA')
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylim([0.8 1.2])
for i=1:length(at_p_val_residue_asa_relative)
    text(i,0.85,num2str(at_p_val_residue_asa_relative(i)))
end



subplot(2,2,4)
hold on
bar(at_mean_residue_neighbor_relative,'BaseValue',1)
title('neighbors')
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylim([0.85 1.15])
for i=1:length(at_p_val_residue_neighbor_relative)
    text(i,0.92,num2str(at_p_val_residue_neighbor_relative(i)))
end




set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;








%gnomad

gnomad_data=readtable([dependency_directory 'human-data/misGnomadScrape.csv']);

hs_af=gnomad_data.Var9;

hs_asa=gnomad_data.sasa;
hs_neighbor=gnomad_data.neighbors;

[~,hs_structure]=structure_types(gnomad_data.secondary);
hs_residue=gnomad_data.refAa;





figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
histogram(hs_af,0:2e3:2.5e5)
set(gca,'YScale','log')
xlabel('gnomad AC')
ylabel('frequency')
axis square



thresh_to_test=10.^(0:0.1:5);

hs_asa_p_val=nan(length(thresh_to_test),1);
hs_asa_effect_size=nan(length(thresh_to_test),1);

hs_neighbor_p_val=nan(length(thresh_to_test),1);
hs_neighbor_effect_size=nan(length(thresh_to_test),1);
for i=1:length(thresh_to_test)
    
    temp_idx=hs_af<=thresh_to_test(i);
    
    if sum(sum(temp_idx))>0
    
        [h p]=ttest2(hs_asa(temp_idx),hs_asa(~temp_idx));

        hs_asa_p_val(i)=p;
        hs_asa_effect_size(i)=mean(hs_asa(temp_idx),'omitnan')./...
            mean(hs_asa(~temp_idx),'omitnan');
        
        
    
        [h p]=ttest2(hs_neighbor(temp_idx),hs_neighbor(~temp_idx));

        hs_neighbor_p_val(i)=p;
        hs_neighbor_effect_size(i)=mean(hs_neighbor(temp_idx),'omitnan')./...
            mean(hs_neighbor(~temp_idx),'omitnan');
        
    end
    
end


subplot(2,4,2)
hold on
scatter(hs_asa_effect_size,-log10(hs_asa_p_val),50,'.k')
scatter(hs_neighbor_effect_size,-log10(hs_neighbor_p_val),25,'ok')
xlim([0.7 1.3])
xlabel('effect size (rare/common)')
ylabel('-log10 p')
ylim([0 350])
legend({'ASA','neighbors'})
axis square
v1=hs_asa_effect_size;
v2=-log10(hs_asa_p_val);
for i=1:length(thresh_to_test)
    text(v1(i),v2(i),num2str(thresh_to_test(i)))
end


%use 200; pretty close to 0.1% here, as well



hs_af_thresh=200;
temp_idx1=hs_af<=hs_af_thresh;

%sum(temp_idx1)
%sum(~temp_idx1)

for i=1:length(structure_labels)
    
    temp_idx2=hs_structure==i;
    
    v1=hs_asa(logical(temp_idx1.*temp_idx2));
    v2=hs_asa(logical(~temp_idx1.*temp_idx2));
    
    hs_mean_structure_asa_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_asa_relative(i)=p;
    
    
    v1=hs_neighbor(logical(temp_idx1.*temp_idx2));
    v2=hs_neighbor(logical(~temp_idx1.*temp_idx2));
    
    hs_mean_structure_neighbor_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_neighbor_relative(i)=p;

end

subplot(2,4,3)
hold on
bar(hs_mean_structure_asa_relative,'BaseValue',1)
title('ASA')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.7 1.2])
for i=1:length(hs_p_val_structure_asa_relative)
    text(i,0.85,num2str(hs_p_val_structure_asa_relative(i)))
end
axis square


subplot(2,4,4)
hold on
bar(hs_mean_structure_neighbor_relative,'BaseValue',1)
title('neighbors')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.9 1.1])
for i=1:length(hs_p_val_structure_neighbor_relative)
    text(i,0.92,num2str(hs_p_val_structure_neighbor_relative(i)))
end
axis square




for i=1:length(aa_labels)
    
    temp_idx2=ismember(hs_residue,aa_labels{i});
    
    v1=hs_asa(logical(temp_idx1.*temp_idx2));
    v2=hs_asa(logical(~temp_idx1.*temp_idx2));
    
    hs_mean_residue_asa_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_residue_asa_relative(i)=p;
    
    
    v1=hs_neighbor(logical(temp_idx1.*temp_idx2));
    v2=hs_neighbor(logical(~temp_idx1.*temp_idx2));
    
    hs_mean_residue_neighbor_relative(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_residue_neighbor_relative(i)=p;

end

subplot(2,2,3)
hold on
bar(hs_mean_residue_asa_relative,'BaseValue',1)
title('ASA')
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylim([0.7 1.3])
for i=1:length(hs_p_val_residue_asa_relative)
    text(i,0.85,num2str(hs_p_val_residue_asa_relative(i)))
end



subplot(2,2,4)
hold on
bar(hs_mean_residue_neighbor_relative,'BaseValue',1)
title('neighbors')
xticks(1:length(aa_labels))
xticklabels(aa_labels)
ylim([0.7 1.3])
for i=1:length(hs_p_val_residue_neighbor_relative)
    text(i,0.92,num2str(hs_p_val_residue_neighbor_relative(i)))
end




set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;




%split by olfactory vs essential? haploinsufficient?
hap_table=readtable([dependency_directory '1-s2.0-S2211124722003175-mmc2.xlsx']);
hap_genes=hap_table.Symbol;


%load pLOF data
lof_table=readtable([dependency_directory 'supplementary_dataset_11_full_constraint_metrics.tsv'],...
    'FileType','text');

idx_to_keep=ismember(lof_table.canonical,'true');
lof_table=lof_table(idx_to_keep,:);

lof_table(isnan(lof_table.oe_lof_upper_rank),:)=[];

[~,sort_idx]=sort(lof_table.oe_lof_upper_rank,'descend');




n_bins=5;

%most to least essential quintiles
quintile_size=floor(length(sort_idx)/n_bins);

for i=1:n_bins
    
    temp_idx=(i-1)*quintile_size+(1:quintile_size);
    gene_sets{i}=lof_table.gene(sort_idx(temp_idx));
    
end



%temp_idx3=ismember(gnomad_data.Var14,gene_sets{5});
temp_idx3=ismember(gnomad_data.Var14,hap_genes);
temp_idx4=ismember(gnomad_data.Var14,gene_sets{1});

% temp_median=median(lof_table.oe_lof_upper_rank);
% temp_idx3=ismember(gnomad_data.Var14,lof_table.gene(lof_table.oe_lof_upper_rank<temp_median));
% temp_idx4=ismember(gnomad_data.Var14,lof_table.gene(lof_table.oe_lof_upper_rank>=temp_median));

%scatter gene-wise stats instead of doing split?
%might just be running out of common variants

for i=1:length(structure_labels)
    
    temp_idx2=hs_structure==i;
    
    v1=hs_asa(logical(temp_idx1.*temp_idx2.*temp_idx3));
    v2=hs_asa(logical(~temp_idx1.*temp_idx2.*temp_idx3));
    
    hs_mean_structure_asa_relative_hap(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_asa_relative_hap(i)=p;
    
    v1=hs_asa(logical(temp_idx1.*temp_idx2.*temp_idx4));
    v2=hs_asa(logical(~temp_idx1.*temp_idx2.*temp_idx4));
    
    hs_mean_structure_asa_relative_other(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_asa_relative_other(i)=p;
    
    
    
    v1=hs_neighbor(logical(temp_idx1.*temp_idx2.*temp_idx3));
    v2=hs_neighbor(logical(~temp_idx1.*temp_idx2.*temp_idx3));
    
    hs_mean_structure_neighbor_relative_hap(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_neighbor_relative_hap(i)=p;
    
    
    v1=hs_neighbor(logical(temp_idx1.*temp_idx2.*temp_idx4));
    v2=hs_neighbor(logical(~temp_idx1.*temp_idx2.*temp_idx4));
    
    hs_mean_structure_neighbor_relative_other(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [h p]=ttest2(v1,v2);
    hs_p_val_structure_neighbor_relative_other(i)=p;

end



figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
hold on
bar([hs_mean_structure_asa_relative_hap' hs_mean_structure_asa_relative_other'],'BaseValue',1)
title('ASA')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.5 1.5])
legend({'more ess.','less ess.'})



subplot(2,2,2)
hold on
bar([hs_mean_structure_neighbor_relative_hap' hs_mean_structure_neighbor_relative_other'],'BaseValue',1)
title('neighbors')
xticks(1:length(structure_labels))
xticklabels(structure_labels)
ylim([0.8 1.2])
legend({'more ess.','less ess.'})





set(gcf,'PaperPositionMode','auto')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'structure_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;







toc




