function [] = plot_evo_emergence_rate(dependency_directory)

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

[n_time_points,n_lineages]=size(condition_slope_mat{i});

for i=1:2

    for j=1:n_time_points
    
        v_temp=condition_slope_mat{i}(j,:);

        n_adapted(i,j)=sum(v_temp>evo_thresh_slope1);

    end

end

hold on

v1=n_adapted(1,:);
v2=n_adapted(2,:);
plot(v1,'k')
plot(v2,'k')
scatter(1:length(v1),v1,25,'k','filled')
scatter(1:length(v2),v2,25,'k','filled')
axis square
xlabel('passage')
ylabel('n adapted lineages')

text(2,70,num2str(n_adapted(1,end)+n_adapted(2,end)))
text(2,65,num2str((n_adapted(1,end)+n_adapted(2,end))/2/8))


end


