function [] = plot_asa_rare_common_structure(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

load([dependency_directory 'structure_mat_1K.mat'])

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

af_thresh=0.01;

temp_idx1=af_mat_1K<=af_thresh;

for i=1:length(structure_labels)

    temp_idx2=structure_mat_1K==i;
    
    v1=reshape(asa_mat_1K(logical(temp_idx1.*temp_idx2)),1,[]);
    v2=reshape(asa_mat_1K(logical(~temp_idx1.*temp_idx2)),1,[]);
    
    to_plot(i)=mean(v1,'omitnan')/mean(v2,'omitnan');

    [h p_val(i)]=ttest2(v1,v2);
    
end



hold on
bar(to_plot,'BaseValue',1)
for i=1:length(p_val)
    text(i,1.05,num2str(p_val(i)),'Rotation',45)
end
ylim([0.85 1.15])
title('ASA (Ang.^2)')
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('rare/common')
axis square

end


