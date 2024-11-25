function [] = plot_mutations_per_gene(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



load([dependency_directory 'mutations_per_gene.mat'])
at_n_muts=at_n_muts(2:end);

hold on
histogram(sc_n_muts,0:25:1000,'Normalization','Probability')
histogram(at_n_muts,0:25:1000,'Normalization','Probability')
histogram(hs_n_muts,0:25:1000,'Normalization','Probability')
axis square
xlabel('mutations per gene')
ylabel('norm. freq.')
legend({'S.c','A.t.','H.s.'})
xlim([0 1e3])
%set(gca,'XScale','log')
ylim([0 0.4])
%medians
text(500,0.3,['S.c. = ' num2str(median(sc_n_muts))])
text(500,0.25,['At.. = ' num2str(median(at_n_muts))])
text(500,0.2,['H.s. = ' num2str(median(hs_n_muts))])


end


