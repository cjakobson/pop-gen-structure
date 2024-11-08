function [] = plot_neighbor_sim_1K_structure(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'neighbor_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
for i=1:length(structure_labels)
    
    temp_idx1=structure_mat_sim==i;
    temp_idx2=structure_mat_1K==i;
    
    v1=reshape(neighbor_mat_1K(temp_idx2),1,[]);
    v2=reshape(neighbor_mat_sim(temp_idx1),1,[]);
    
    to_plot(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [~,p_val(i)]=ttest2(v1,v2);
    
end


hold on
bar(to_plot,'BaseValue',1)
for i=1:length(p_val)
    text(i,1.1,num2str(p_val(i)),'Rotation',45)
end
ylim([0.7 1.3])
title('C_\alpha within 10 Ang.')
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('1K normalized to simulated')
axis square


end


