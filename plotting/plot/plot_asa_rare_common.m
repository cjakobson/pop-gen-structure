function [] = plot_asa_rare_common(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])
load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

af_thresh=0.01;

v1=reshape(asa_mat_1K(af_mat_1K<=af_thresh),1,[]);
v1=v1(~isnan(v1));

v2=reshape(asa_mat_1K(af_mat_1K>af_thresh),1,[]);
v2=v2(~isnan(v2));



hold on
histogram(v1,0:10:250,'Normalization','probability')
histogram(v2,0:10:250,'Normalization','probability')
axis square
legend({'rare','common'})
%title('simulated SNPs')
ylabel('relative freq.')
xlabel('ASA (Ang.^2)')
%set(gca,'YScale','log')
[h p]=kstest2(v1,v2);
text(200,0.1,['p = ' num2str(p)])
text(200,0.08,['rare = ' num2str(length(v1))])
text(200,0.06,['common = ' num2str(length(v2))])


end


