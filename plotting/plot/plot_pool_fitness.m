function [] = plot_pool_fitness(dependency_directory)

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


%start_time=24;

%evo_thresh_slope1=0.75;
%evo_thresh_slope2=0.5;


v_row={'A','B','C','D','E','F','G','H',...
    'I','J','K','L','M','N','O','P'};
%get pool indices wrt evo plates
pool_data=readtable([dependency_directory '20231211_fpr1_pools.xlsx']);
pools_to_use=[1:3 7:9 13:15 19:21];

pool_plate_idx=pool_data.x96wellPlate(pools_to_use);

pool_temp_start=pool_data.x96wellStart(pools_to_use);
pool_temp_end=pool_data.x96wellEnd(pools_to_use);

for i=1:length(pool_temp_start)
    
    pool_row{i}=pool_temp_start{i}(1);
    
    pool_column_start(i)=str2num(pool_temp_start{i}(2:end));
    pool_column_end(i)=str2num(pool_temp_end{i}(2:end));
    
end

%look up evo plate positions
pool_coordinates=readtable([dependency_directory '20230427_rearray_96well_coordinates.csv']);

for i=1:length(pool_plate_idx)
    
    plate_idx=pool_coordinates.Var8==pool_plate_idx(i);
    
    row_idx=ismember(pool_coordinates.Var9,pool_row{i});
    
    start_idx=pool_coordinates.Var10>=pool_column_start(i);
    end_idx=pool_coordinates.Var10<=pool_column_end(i);
    
    pool_coord_idx{i}=find(plate_idx.*row_idx.*start_idx.*end_idx);
    
end

v_row={'A','B','C','D','E','F','G','H',...
    'I','J','K','L','M','N','O','P',...
    'Q','R','S','T','U','V','W','X',...
    'Y','Z','AA','AB','AC','AD','AE','AF'};
%retrieve final fitness
adapted_thresh=0.75;
for i=1:length(pool_coord_idx)
    
    temp_coord_idx=pool_coord_idx{i};
    
    for j=1:length(temp_coord_idx)
        
        temp_evo_plate=str2num(pool_coordinates.Var1{temp_coord_idx(j)}(end));
        
        temp_evo_row=find(ismember(v_row,pool_coordinates.Var2{temp_coord_idx(j)}));
        
        temp_evo_col=pool_coordinates.Var3(temp_coord_idx(j));
        
        temp_evo_absolute_index=48*(temp_evo_row-1)+temp_evo_col;
        
        pool_absolute_idx{i}(j)=temp_evo_absolute_index;
        
        %v_final=condition_slope_mat{temp_evo_plate}(end,:);
        v_fitness=phenotyping_mean_rap(temp_evo_plate,:);
        pool_fitness{i}(j)=v_fitness(pool_absolute_idx{i}(j));
        
    end
    
    n_adapted(i)=sum(pool_fitness{i}>adapted_thresh);
    
end

hold on
easy_box_with_dots(pool_fitness)
ylim([0 2])
ylabel('growth rate')
plot(xlim,[adapted_thresh adapted_thresh],':r')
axis square

%sum(cellfun(@length,pool_coord_idx))
%sum(n_adapted)

end


