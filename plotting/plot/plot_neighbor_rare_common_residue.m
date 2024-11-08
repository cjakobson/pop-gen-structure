function [] = plot_neighbor_rare_common_residue(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

load([dependency_directory 'residue_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


af_thresh=0.01;

temp_idx1=af_mat_1K<=af_thresh;

for i=1:length(aa_labels)

    temp_idx2=residue_mat_1K==i;
    
    v1=reshape(neighbor_mat_1K(logical(temp_idx1.*temp_idx2)),1,[]);
    v2=reshape(neighbor_mat_1K(logical(~temp_idx1.*temp_idx2)),1,[]);
    
    to_plot(i)=mean(v1,'omitnan')/mean(v2,'omitnan');

    [h p_val(i)]=ttest2(v1,v2);

    
end



hold on
bar(to_plot,'BaseValue',1)
for i=1:length(p_val)
    text(i,0.98,num2str(p_val(i)),'Rotation',-45)
end
ylim([0.85 1.15])
title('C_\alpha within 10 Ang.')
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
ylabel('rare/common')
axis square
xlim([0.5 length(aa_labels)+0.5])

end


