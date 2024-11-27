function [] = plot_evo_rates_rap(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



exp_name={'round1','round2','round3','round4'};

plate_labels={'rap1','rap2','YPDonly'};

condition_names={'rap','YPD'};


for i=1:length(exp_name)
    
    for j=1:length(plate_labels)
        
        if i<=2
            k=1;
        else
            k=2;
        end
    

        phenotyping_slope_mat{i,j}=table2array(readtable([dependency_directory 'gitter-data/'...
            exp_name{i} '_' plate_labels{j} '=>' condition_names{k} '_phen_slope.csv']));

        phenotyping_slope_names{i,j}=[exp_name{i} '_' plate_labels{j} '=>' condition_names{k}];

    end

end



%scatter replicates
v3=[];
v4=[];
for i=1:2%length(plate_labels)

    v1=phenotyping_slope_mat{1,i};
    v2=phenotyping_slope_mat{2,i};

    phenotyping_mean_rap(i,:)=mean([v1; v2]);

    

end

for i=1:2%length(plate_labels)

    v1=phenotyping_slope_mat{3,i};
    v2=phenotyping_slope_mat{4,i};

    phenotyping_mean_ypd(i,:)=mean([v1; v2]);

end





adapted_slope_thresh=0.75;


temp_labels={'in rap','in rap','in YPD','in YPD'};

for i=1:4
    to_plot{i}=[];
end

for i=1:2%length(plate_labels)

    v1=phenotyping_mean_rap(i,:);
    v2=phenotyping_mean_ypd(i,:);

    adapted_idx=v1>adapted_slope_thresh;

    %clear to_plot

    to_plot{1}=[to_plot{1}; v1(~adapted_idx)'];
    to_plot{2}=[to_plot{2}; v1(adapted_idx)'];

    to_plot{3}=[to_plot{3}; v2(~adapted_idx)'];
    to_plot{4}=[to_plot{4}; v2(adapted_idx)'];

end


easy_box(to_plot)
title(plate_labels{i})

xticks(1:length(to_plot))
xtickangle(45)
xticklabels(temp_labels)
ylabel('growth rate')
ylim([0 2])
axis square
for j=1:length(to_plot)
    text(j,0.1,num2str(length(to_plot{j})))
end
for i=1:2
    [h p]=ttest2(to_plot{2*(i-1)+1},to_plot{2*(i-1)+2});
    text(2*(i-1)+1.5,1.8,num2str(p))
end


end


