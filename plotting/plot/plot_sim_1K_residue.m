function [] = plot_sim_1K_residue(dependency_directory)

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
    
    v1=reshape(asa_mat_1K(temp_idx2),1,[]);
    v2=reshape(asa_mat_sim(temp_idx1),1,[]);
    
    to_plot1(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [~,p_val1(i)]=ttest2(v1,v2);
    
end



load([dependency_directory 'neighbor_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])


for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    v1=reshape(neighbor_mat_1K(temp_idx2),1,[]);
    v2=reshape(neighbor_mat_sim(temp_idx1),1,[]);
    
    to_plot2(i)=mean(v1,'omitnan')/mean(v2,'omitnan');
    
    [~,p_val2(i)]=ttest2(v1,v2);
    
end





hold on
%bar(to_plot,'BaseValue',1)
bar([to_plot1; to_plot2]','BaseValue',1)
for i=1:length(p_val1)
    text(i,to_plot1(i),num2str(p_val1(i)),'Rotation',45)
    text(i,to_plot2(i),num2str(p_val2(i)),'Rotation',-45)
end
ylim([0.6 1.8])
%title('ASA (Ang.^2)')
legend({'ASA (Ang.^2)','C_\alpha within 10 Ang.'})
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
xlim([0 length(aa_labels)+1])
ylabel('1K normalized to simulated')
%axis square


end


