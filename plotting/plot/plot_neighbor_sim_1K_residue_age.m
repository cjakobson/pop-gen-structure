function [] = plot_neighbor_sim_1K_residue_age(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'neighbor_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory 'residue_mat_sim.mat'])
load([dependency_directory 'residue_mat_1K.mat'])

%filter on genes with too few 1K variants
cov_thresh=50;
low_coverage_idx=sum(~isnan(neighbor_mat_1K),2)<cov_thresh;

neighbor_mat_1K(low_coverage_idx,:)=nan;

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

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

%young_idx=gene_age==6;
young_idx=gene_age>=2;
%old_idx=gene_age==1;
old_idx=gene_age<=1;




gene_residue_neighbor_mat=nan(length(genes_to_use),length(aa_labels));
for i=1:length(aa_labels)
    
    for j=1:length(genes_to_use)
        
        temp_structure_idx_1K=residue_mat_1K(j,:)==i;
        temp_structure_idx_sim=residue_mat_sim(j,:)==i;
        
        temp_neighbor_1K=neighbor_mat_1K(j,temp_structure_idx_1K);
        temp_neighbor_sim=neighbor_mat_sim(j,temp_structure_idx_sim);
        
        gene_residue_neighbor_mat(j,i)=mean(temp_neighbor_1K,'omitnan')/...
            mean(temp_neighbor_sim,'omitnan');
    
    end
    
end


for i=1:length(aa_labels)
    
    v1=gene_residue_neighbor_mat(old_idx,i);
    v2=gene_residue_neighbor_mat(young_idx,i);

    to_plot1(i)=mean(v1,'omitnan');
    to_plot2(i)=mean(v2,'omitnan');
    
    [h,p_val(i)]=ttest2(v1,v2);
    
end


hold on
bar([to_plot1; to_plot2]','BaseValue',1)
for i=1:length(p_val)
    text(i,1.1,num2str(p_val(i)),'Rotation',45)
end
ylim([0.8 1.15])
xlim([0.5 length(to_plot1)+0.5])
title('C_\alpha within 10 Ang.')
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
ylabel('1K normalized to simulated')
legend({'ancient','youngest'})
axis square


end


