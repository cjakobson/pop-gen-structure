function [] = plot_wgs_clones(dependency_directory,evo_plate_pos)



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

v_row={'A','B','C','D','E','F','G','H',...
    'I','J','K','L','M','N','O','P'};
for i=1:length(evo_plate_pos)
    
    plate_idx(i)=str2num(evo_plate_pos{i}(1));
    
    temp_row=evo_plate_pos{i}(2);
    row_idx(i)=find(ismember(v_row,temp_row));
    
    col_idx(i)=str2num(evo_plate_pos{i}(3:end));
    
    absolute_idx(i)=48*(row_idx(i)-1)+col_idx(i);
    
end

hold on
for i=1:length(evo_plate_pos)
    
    [n_time_points,~]=size(condition_spot_size_mat{plate_idx(i)});

    v1=1:n_time_points;
    v2=condition_spot_size_mat{plate_idx(i)}(:,absolute_idx(i));

    plot(v1(start_time:end),v2(start_time:end),'-r')
    scatter(v1(start_time:end),v2(start_time:end),20,'r','filled')
    
end



ylim([0 3000])
%title(condition_names{i})
ylabel('spot size')
xlabel('time (hrs)')





end


