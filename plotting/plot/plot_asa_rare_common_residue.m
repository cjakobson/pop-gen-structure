function [] = plot_asa_rare_common_residue(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

load([dependency_directory 'residue_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


af_thresh=0.01;

temp_idx1=af_mat_1K<=af_thresh;

for i=1:length(aa_labels)

    temp_idx2=residue_mat_1K==i;
    
    to_plot(i)=mean(mean(asa_mat_1K(logical(temp_idx1.*temp_idx2)),'omitnan'),'omitnan')/...
        mean(mean(asa_mat_1K(logical(~temp_idx1.*temp_idx2)),'omitnan'),'omitnan');

    
end



hold on
bar(to_plot,'BaseValue',1)
ylim([0.7 1.2])
title('ASA (Ang.^2)')
xticks(1:length(aa_labels))
xtickangle(45)
xticklabels(aa_labels)
ylabel('rare/common')
axis square
xlim([0.5 length(aa_labels)+0.5])

end


