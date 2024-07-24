
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];
output_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-output/'];

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

load([dependency_directory 'simulated_mutation_positions.mat'])

load([dependency_directory '1K_mutation_positions.mat'])
load([dependency_directory '1K_mutation_af.mat'])


structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','R','H','K','D','E','C','G','P'};



%general histograms/correlations for structure stats (just for all
%residues)




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
ylim([0.6 1.4])
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
print([output_directory 'figure_1'],'-dsvg','-r0')
print([output_directory 'figure_1'],'-djpeg','-r300')







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
    
    asa_1K(i)=mean(temp_asa(pos_1K));
    asa_sim(i)=mean(temp_asa(pos_sim));
    rel_asa(i)=asa_1K(i)/asa_sim(i);
    
    neighbor_1K(i)=mean(temp_neighbor(pos_1K));
    neighbor_sim(i)=mean(temp_neighbor(pos_sim));
    rel_neighbor(i)=neighbor_1K(i)/neighbor_sim(i);
    
    %also how many mutations fixed relative to all possible
    %to control for general essentiality
    mutations_1K(i)=length(pos_1K)/protein_length(i);
    mutations_sim(i)=length(pos_sim)/protein_length(i);
    rel_mutations(i)=length(pos_1K)/length(pos_sim);
    
end



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
hold on
histogram(asa_1K)%,'Normalization','probability')
histogram(asa_sim)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('ASA')
ylabel('number of genes')



subplot(2,4,2)
hold on
histogram(neighbor_1K)%,'Normalization','probability')
histogram(neighbor_sim)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('neighbors')
ylabel('number of genes')


subplot(2,4,3)
hold on
histogram(mutations_1K)%,'Normalization','probability')
histogram(mutations_sim)%,'Normalization','probability')
axis square
legend({'1K','sim'})
title('mutations/residue')
ylabel('number of genes')





subplot(2,4,5)
hold on
scatter(rel_asa,rel_neighbor,10,'k','filled')
axis square
xlabel('relative ASA 1K/sim')
ylabel('relative neighbors 1K/sim')


subplot(2,4,6)
hold on
scatter(rel_mutations,rel_asa,10,'k','filled')
axis square
xlabel('relative # mutations 1K/sim')
ylabel('relative ASA 1K/sim')


subplot(2,4,7)
hold on
scatter(rel_mutations,rel_neighbor,10,'k','filled')
axis square
xlabel('relative # mutations 1K/sim')
ylabel('relative neighbors 1K/sim')




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2'],'-dsvg','-r0')
print([output_directory 'figure_2'],'-djpeg','-r300')







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




