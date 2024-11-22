function [] = plot_niche_example(dependency_directory,chr_to_plot,pos_to_plot)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

niche_input=readtable([dependency_directory '1K_common_annotated_niche.csv']);

niche_names={'Baking','Dairy','Ferm','Human','Other','Plant','Soil','Waste'};

idx_to_plot=find(ismember(niche_input.chr,chr_to_plot).*niche_input.pos==pos_to_plot);

v1=table2array(niche_input(idx_to_plot,12:19)); %ref
v2=table2array(niche_input(idx_to_plot,20:27)); %   alt

v3=v2./(v1+v2);

hold on
scatter(1:length(v3),v3,150,'k','filled')
ylim([0 1])
xlim([0.5 length(niche_names)+0.5])
xticks(1:length(niche_names))
xtickangle(45)
xticklabels(niche_names)
title(niche_input.gene{idx_to_plot})
ylabel('alt allele freq.')
text(niche_input.niche_id(idx_to_plot),1,num2str(niche_input.niche_q_value(idx_to_plot)))

end


