function [] = plot_asa_sim_1K_residue(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'residue_mat_sim.mat'])
load([dependency_directory 'residue_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

 
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    to_plot(i)=mean(mean(asa_mat_1K(temp_idx2),'omitnan'),'omitnan')/...
        mean(mean(asa_mat_sim(temp_idx1),'omitnan'),'omitnan');
    
end


hold on
bar(to_plot,'BaseValue',1)
ylim([0.7 1.8])
title('ASA (Ang.^2)')
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
xlim([0 length(aa_labels)+1])
ylabel('1K normalized to simulated')
axis square


end


