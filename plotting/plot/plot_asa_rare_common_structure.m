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
    
    to_plot(i)=mean(mean(asa_mat_1K(logical(temp_idx1.*temp_idx2)),'omitnan'),'omitnan')/...
        mean(mean(asa_mat_1K(logical(~temp_idx1.*temp_idx2)),'omitnan'),'omitnan');

    
end



hold on
bar(to_plot,'BaseValue',1)
ylim([0.85 1.15])
title('ASA (Ang.^2)')
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('rare/common')
axis square

end


