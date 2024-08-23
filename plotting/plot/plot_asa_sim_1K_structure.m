function [] = plot_asa_sim_1K_structure(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
for i=1:length(structure_labels)
    
    temp_idx1=structure_mat_sim==i;
    temp_idx2=structure_mat_1K==i;
    
    to_plot(i)=mean(mean(asa_mat_1K(temp_idx2),'omitnan'),'omitnan')/...
        mean(mean(asa_mat_sim(temp_idx1),'omitnan'),'omitnan');
    
end


hold on
bar(to_plot,'BaseValue',1)
ylim([0.7 1.3])
title('ASA (Ang.^2)')
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('1K normalized to simulated')
axis square


end


