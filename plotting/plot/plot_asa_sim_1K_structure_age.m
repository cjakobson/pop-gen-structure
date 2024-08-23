function [] = plot_neighbor_sim_1K_structure_age(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
load([dependency_directory 'gene_names.mat'])
 
 
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




gene_structure_asa_mat=nan(length(genes_to_use),length(structure_labels));
for i=1:length(structure_labels)
    
    for j=1:length(genes_to_use)
        
        temp_structure_idx_1K=structure_mat_1K(j,:)==i;
        temp_structure_idx_sim=structure_mat_sim(j,:)==i;
        
        temp_asa_1K=asa_mat_1K(j,temp_structure_idx_1K);
        temp_asa_sim=asa_mat_sim(j,temp_structure_idx_sim);
        
        gene_structure_asa_mat(j,i)=mean(temp_asa_1K,'omitnan')/...
            mean(temp_asa_sim,'omitnan');
    
    end
    
end

to_plot1=mean(gene_structure_asa_mat(old_idx,:),'omitnan');
to_plot2=mean(gene_structure_asa_mat(young_idx,:),'omitnan');




hold on
bar([to_plot1; to_plot2]','BaseValue',1)
ylim([0.8 1.4])
title('ASA (Ang.^2) 1K/sim')
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('1K normalized to simulated')
legend({'ancient','youngest'})
axis square


end


