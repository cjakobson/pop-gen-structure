function [] = plot_neighbor_sim_raw_residue(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'neighbor_mat_sim.mat'])

load([dependency_directory 'residue_mat_sim.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

 
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    
    to_plot(i)=mean(mean(neighbor_mat_sim(temp_idx1),'omitnan'),'omitnan');
    
end


hold on
bar(to_plot)
ylim([0 20])
title('C_\alpha within 10 Ang.')
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
xlim([0 length(aa_labels)+1])
ylabel('simulated')
axis square


end


