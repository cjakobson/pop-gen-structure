function [] = plot_structure_all(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

load([dependency_directory 'structure_mat_1K.mat'])

for i=1:length(structure_labels)

    v1(i)=sum(sum(structure_mat_1K==i,'omitnan'),'omitnan');

end

n_total=sum(sum(~isnan(structure_mat_1K),'omitnan'),'omitnan');
v1=v1./n_total;




at_data=readtable([dependency_directory 'arabidopsis-data/1001misAnnotated.csv']);
v2=structure_types(at_data.secondary);


hs_data=readtable([dependency_directory 'human-data/misGnomadScrape.csv']);
v3=structure_types(hs_data.secondary);


hold on
bar([v1;v2;v3]')
ylim([0 0.5])
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('norm. frequency')
axis square
legend({'S.c','A.t.','H.s.'})


end


