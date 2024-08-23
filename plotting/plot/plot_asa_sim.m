function [] = plot_asa_sim(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_data.mat'])

v1=reshape(asa_mat,1,[]);
v1=v1(~isnan(v1));


hold on
histogram(v1,0:10:250)%,'Normalization','probability')
axis square
legend({'all residues'})
title('all residues in proteome')
ylabel('number of residues')
xlabel('ASA (Ang.^2)')
set(gca,'YScale','log')

end


