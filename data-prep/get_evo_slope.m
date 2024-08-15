
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])



exp_name='evolution';

plate_names={'a','b','c','d','e'};

n_scanners=1;

condition_names={'YPD+rap1','YPD+rap2','YPDonly'};

%starting date and time of the exp
start_date=20230306;

start_time=9*60+52+30/60;   %in mins



initial_offset=6; %wait 6hrs after pinning
total_time=12;   %measure rate over next 12h interval
total_offset=24; %remeasure every 24h

n_passages=9;




for i=1:length(condition_names)
    
    i
    
    input_table=readtable([dependency_directory 'gitter-data/'...
        condition_names{i} '_evo.csv']);
    
    growth_mat=table2array(input_table(:,2:end));
    v_time=table2array(input_table(:,1));
    
    [n_time_points,n_strains]=size(growth_mat);
    
    for j=1:n_passages
                
        j
        
        for k=1:n_strains

            min_time=(initial_offset+(total_offset*(j-1)))*60;
            max_time=(initial_offset+(total_offset*(j-1))+total_time)*60;

            temp_idx=logical((v_time>min_time).*(v_time<max_time));

            temp_lm=fitlm(v_time(temp_idx),growth_mat(temp_idx,k));

            slope_mat{i}(j,k)=table2array(temp_lm.Coefficients(2,1));

        end
            
    end
    
    to_output=table(slope_mat{i});
    writetable(to_output,[dependency_directory 'gitter-data/'...
        condition_names{i} '_evo_slope.csv']);
    
end



