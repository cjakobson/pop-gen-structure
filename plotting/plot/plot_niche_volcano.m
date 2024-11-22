function [] = plot_niche_volcano(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

niche_input=readtable([dependency_directory '1K_common_annotated_niche.csv']);

v1=log10(niche_input.niche_enrichment);
v2=niche_input.niche_q_value;

hold on
scatter(v1,v2,10,'k','filled')
axis square
xlabel('log_{10} alt allele enrichment')
ylabel('log_{10}q')
xlim([-3 3])
ylim([-5 25])
plot(xlim,[2 2],':r')
text(2,-3,['n sig = ' num2str(sum(v2>2))])
text(2,-4,['n total = ' num2str(length(v2))])


end


