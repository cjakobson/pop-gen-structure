function [] = plot_neighbor_rare_common(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])
load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

af_thresh=0.01;

v1=reshape(neighbor_mat_1K(af_mat_1K<=af_thresh),1,[]);
v1=v1(~isnan(v1));

v2=reshape(neighbor_mat_1K(af_mat_1K>af_thresh),1,[]);
v2=v2(~isnan(v2));



hold on
histogram(v1,0:2:50,'Normalization','probability')
histogram(v2,0:2:50,'Normalization','probability')
axis square
legend({'rare','common'})
%title('simulated SNPs')
ylabel('relative freq.')
xlabel('C_\alpha within 10 Ang.')
%set(gca,'YScale','log')

end


