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

    
    hold on
    scatter(v1,v2,10,'k','filled')
    %title([plate_labels{i} ' in ' condition_names{1}])
    xlim([0 2])
    ylim(xlim)
    axis square
    xlabel('replicate 1')
    ylabel('replicate 2')

    v3=[v3;v1'];
    v4=[v4;v2'];

end

[r p]=corr(v3,v4,'rows','complete');
text(1,0.5,['r = ' num2str(r)])

adapted_slope_thresh=0.75;
plot([adapted_slope_thresh adapted_slope_thresh],ylim,':k')
plot(xlim,[adapted_slope_thresh adapted_slope_thresh],':k')

end


