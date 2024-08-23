function [] = scatter_asa_neighbor_gene(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'neighbor_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory 'protein_length.mat'])
load([dependency_directory 'gene_names.mat'])

for i=1:length(genes_to_use)
    
    temp_asa_sim=asa_mat_sim(i,:);
    temp_asa_1K=asa_mat_1K(i,:);
    
    v1(i)=mean(temp_asa_1K,'omitnan')/mean(temp_asa_sim,'omitnan');
    
    
    temp_neighbor_sim=neighbor_mat_sim(i,:);
    temp_neighbor_1K=neighbor_mat_1K(i,:);
    
    v2(i)=mean(temp_neighbor_1K,'omitnan')/mean(temp_neighbor_sim,'omitnan');
    
    
    temp_n_mutations_sim=sum(~isnan(temp_neighbor_sim));
    temp_n_mutations_1K=sum(~isnan(temp_neighbor_1K));
    
    v3(i)=mean(temp_n_mutations_1K,'omitnan')/mean(temp_n_mutations_sim,'omitnan');
    
end


hold on
scatter(v1,v2,5,'k','filled')
axis square
xlabel('ASA (Ang.^2) 1K/sim')
ylabel('C_\alpha within 10 Ang. 1K/sim')
xlim([0.5 2])
ylim([0.5 1.2])
[r p]=corr(v1',v2','rows','complete');
text(1.5,1.15,num2str(r))
text(1.5,1.1,num2str(p))
% bar(to_plot,'BaseValue',1)
% ylim([0.7 1.3])
% title('ASA (Ang.^2)')
% xticks(1:length(structure_labels))
% xtickangle(45)
% xticklabels(structure_labels)
% ylabel('1K normalized to simulated')
% axis square


end


