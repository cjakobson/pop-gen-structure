
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

plate_labels={'rap1','rap2','NA','NA','YPDonly'};

for i=1:length(exp_name)
    
    
    if i<=2
        condition_names={'rap1=>rap','rap2=>rap',...
            'NA','NA','YPDonly=>rap',...
            'NA','NA','NA','NA','NA',};
    else
        condition_names={'rap1=>YPD','rap2=>YPD',...
            'NA','NA',...
            'YPDonly=>YPD',...
            'NA','NA','NA','NA','NA',};
    end
    
    %get filenames
    dir_base=[dependency_directory 'gitter-data/'];

    to_read=dir([dir_base 'phenotyping/' exp_name{i}]);
    
    m=1;
    
    for j=1:length(to_read)
        j
        temp_name=to_read(j).name;

        if length(temp_name)>3

            if strcmp(temp_name((end-2):end),'dat')
                
                %grab dates and times
                temp_str=strsplit(temp_name,'_');
                dates{i}(m)=str2num(temp_str{1});
                times{i}{m}=temp_str{end-2};

                %convert time to minutes and normalize dates
                temp_date=dates{i}(m)-start_date(i);

                temp_hours=times{i}{m}(1:2);
                temp_mins=times{i}{m}(3:4);
                temp_secs=times{i}{m}(5:6);


                time_mins{i}(m)=str2num(temp_hours)*60+str2num(temp_mins)+...
                    str2num(temp_secs)/60+24*60*temp_date-start_time(i);

                temp_str2=strsplit(temp_str{end},'.');

                scanner{i}(m)=str2num(temp_str2{1}(1));
                plates{i}{m}=temp_str2{1}(2:end);

                %get sga data
                sga_mat{i}{m}=readtable([dir_base 'phenotyping/' exp_name{i} '/' temp_name]);

                mat_to_process{i}(:,m)=sga_mat{i}{m}.size;

                m=m+1;

                
            end
            
        end
        
    end
    
    
    plates_to_output={'a','b','e'};


    for j=1:length(plates_to_output)

        idx_to_use=ismember(plates{i},plates_to_output{j});

        mat_to_output=mat_to_process{i}(:,idx_to_use);
        v_time=time_mins{i}(idx_to_use);

        to_output=table(v_time',mat_to_output');
        writetable(to_output,[dependency_directory 'gitter-data/'...
            exp_name{i} '_' condition_names{ismember(plate_names,plates_to_output{j})} '_phen.csv'])

    end

end






