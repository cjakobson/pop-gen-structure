function [] = plot_evo_example_curves(dependency_directory)

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

evo_thresh_slope1=1;
evo_thresh_slope2=0.5;


for i=1:2%length(condition_names)

    [n_time_points,~]=size(condition_spot_size_mat{i});

    v_final=condition_slope_mat{i}(end,:);
    adapted_idx{i}=find(v_final>evo_thresh_slope1);
    non_adapted_idx{i}=find(v_final<=evo_thresh_slope2);
    
    hold on

    v1=1:n_time_points;
    v2=condition_spot_size_mat{i}(:,non_adapted_idx{i}(1));

    plot(v1(start_time:end),v2(start_time:end),'-k')
    scatter(v1(start_time:end),v2(start_time:end),20,'k','filled')



    v1=1:n_time_points;
    v2=condition_spot_size_mat{i}(:,adapted_idx{i}(11));

    plot(v1(start_time:end),v2(start_time:end),'-r')
    scatter(v1(start_time:end),v2(start_time:end),20,'r','filled')
    
    ylim([0 3000])
    %title(condition_names{i})
    ylabel('spot size')
    xlabel('time (hrs)')


end





end


