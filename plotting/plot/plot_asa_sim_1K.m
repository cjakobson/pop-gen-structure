function [] = plot_asa_sim_1K(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

v1=reshape(asa_mat_sim,1,[]);
v1=v1(~isnan(v1));
v2=reshape(asa_mat_1K,1,[]);
v2=v2(~isnan(v2));


hold on
histogram(v1,0:10:250,'Normalization','probability')
histogram(v2,0:10:250,'Normalization','probability')
axis square
legend({'simulated SNPs','1K SNPs'})
%title('simulated SNPs')
ylabel('number of residues')
xlabel('ASA (Ang.^2)')
%set(gca,'YScale','log')

end

