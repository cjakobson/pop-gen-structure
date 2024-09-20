function [] = plot_neighbor_sim(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'neighbor_data.mat'])

v1=reshape(neighbor_mat,1,[]);
v1=v1(~isnan(v1));


hold on
histogram(v1,0:2:50)%,'Normalization','probability')
axis square
legend({'all residues'})
title('all residues in proteome')
ylabel('number of residues')
xlabel('C_\alpha within 10 Ang.')
set(gca,'YScale','log')
ylim([1e3 1e6])

end


