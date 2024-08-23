function [] = plot_asa_sim_1K_superfam(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'domain_mat_sim.mat'])
load([dependency_directory 'domain_mat_1K.mat'])

to_plot{1}=asa_mat_1K(domain_mat_1K==1);
to_plot{3}=asa_mat_1K(domain_mat_1K==0);

to_plot{2}=asa_mat_sim(domain_mat_sim==1);
to_plot{4}=asa_mat_sim(domain_mat_sim==0);


hold on
easy_box(to_plot)
axis square
temp_labels={'1K in domain','sim in domain',...
    '1K outside domain','sim outside domain'};
xticks(1:length(temp_labels))
xtickangle(45)
xticklabels(temp_labels)
ylabel('ASA (Ang.^2)')
ylim([0 300])

end


