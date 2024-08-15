
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])



exp_name={'round1','round2','round3','round4'};
start_date=[20230323,20230324,20230325,20230326];
start_time=[11*60+53+13/60,10*60+34+54/60,11*60+54+23/60,10*60+27+53/60];

plate_names={'a','b','c','d','e'};

n_scanners=2;

plate_labels={'rap1','rap2','YPDonly'};

condition_names={'rap','YPD'};



initial_offset=6; %wait 6hrs after pinning
total_time=12;   %measure rate over next 12h interval
total_offset=24; %remeasure every 24h

n_passages=9;




for i=1:length(exp_name)
    
    for j=1:length(plate_labels)
        
        if i<=2
            k=1;
        else
            k=2;
        end
    
        input_table=readtable([dependency_directory 'gitter-data/'...
            exp_name{i} '_' plate_labels{j} '=>' condition_names{k} '_phen.csv']);


        growth_mat=table2array(input_table(:,2:end));
        v_time=table2array(input_table(:,1));

        [n_time_points,n_strains]=size(growth_mat);

        for l=1:n_strains

            min_time=(initial_offset)*60;
            max_time=(initial_offset+total_time)*60;

            temp_idx=logical((v_time>min_time).*(v_time<max_time));

            temp_lm=fitlm(v_time(temp_idx),growth_mat(temp_idx,l));

            slope_mat{i,j,k}(l)=table2array(temp_lm.Coefficients(2,1));

        end

        to_output=table(slope_mat{i,j,k});
        writetable(to_output,[dependency_directory 'gitter-data/'...
            exp_name{i} '_' plate_labels{j} '=>' condition_names{k} '_phen_slope.csv'])


        
    end
    
end
    
    
    




