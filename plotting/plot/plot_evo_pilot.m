function [] = plot_evo_pilot(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



%dotplot of final adapted growth rates from pilot
strains={'BY4741','BY4741Dmsh6','BY4743'};%,'BY4743DDmsh6'};

clear to_plot
for i=1:length(strains)

    spot_size_table{i}=readtable([dependency_directory 'gitter-data/' strains{i} '_pilot.csv']);
    spot_size_mat{i}=table2array(spot_size_table{i}(:,2:end));

    to_plot{i}=spot_size_mat{i}(end,:);

end





hold on
easy_dotplot(to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(strains)
ylabel('final spot size')
%axis square
title('rapamycin pilot')
xlim([0.5 3.5])

evo_thresh_spot_size=1000;%1500;

plot(xlim,[evo_thresh_spot_size evo_thresh_spot_size],':r')
for i=1:length(to_plot)
    n_adapted=sum(to_plot{i}>evo_thresh_spot_size);
    text(i,2500,num2str(n_adapted),'Rotation',45)
    f_adapted=sum(to_plot{i}>evo_thresh_spot_size)/length(to_plot{i});
    text(i,2750,num2str(f_adapted),'Rotation',45)
end



end


