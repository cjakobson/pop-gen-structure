function [] = plot_evo_rates(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%plot example curves over time
condition_names={'YPD+rap1','YPD+rap2','YPDonly'};
start_passage=2;

for i=1:length(condition_names)
    
    condition_spot_size_table{i}=readtable([dependency_directory 'gitter-data/' condition_names{i} '_evo.csv']);
    condition_spot_size_mat{i}=table2array(condition_spot_size_table{i}(:,2:end));



    condition_slope_table{i}=readtable([dependency_directory 'gitter-data/' condition_names{i} '_evo_slope.csv']);
    condition_slope_mat{i}=table2array(condition_slope_table{i});

    condition_slope_mat{i}=condition_slope_mat{i}(start_passage:end,:);


end




start_time=24;

evo_thresh_slope1=0.75;
evo_thresh_slope2=0.5;


for i=1:2%length(condition_names)
    
    

    [n_time_points,n_strains]=size(condition_slope_mat{i});

    hold on

    v_final=condition_slope_mat{i}(end,:);
    adapted_idx{i}=find(v_final>evo_thresh_slope1);
    non_adapted_idx{i}=find(v_final<=evo_thresh_slope2);

    %example adapted and non-adapted trajectories

    for j=1:length(adapted_idx{i})

        v1=1:n_time_points;
        v2=condition_slope_mat{i}(:,adapted_idx{i}(j));

        plot(v1,v2,'-r')
        scatter(v1,v2,20,'r','filled')

    end

    for j=101:200%length(adapted_idx)

        v1=1:n_time_points;
        v2=condition_slope_mat{i}(:,non_adapted_idx{i}(j));

        plot(v1,v2,'-k')
        scatter(v1,v2,20,'k','filled')

    end

    ylim([0 2.5])
    xlim([0.5 n_time_points+0.5])
    %axis square
    title(condition_names{i})

end
ylabel('growth rate')
xlabel('passage #')

text(1,2.5,['n adapted = ' num2str(length(adapted_idx{1})) ' '...
    num2str(length(adapted_idx{2}))])




end


