
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])



exp_name='evolution';

plate_names={'a','b','c','d','e'};

n_scanners=1;

condition_names={'YPD+rap1','YPD+rap2','NA','NA',...
    'YPDonly'};

%get filenames
dir_base=[dependency_directory 'gitter-data/'];

to_read=dir([dir_base exp_name]);


%starting date and time of the exp
start_date=20230306;

start_time=9*60+52+30/60;   %in mins


%get sga data
m=1;

for i=1:length(to_read)
    i
    temp_name=to_read(i).name;
    
    if length(temp_name)>3
        
        if strcmp(temp_name((end-2):end),'dat')
            
            %grab dates and times
            temp_str=strsplit(temp_name,'_');
            dates(m)=str2num(temp_str{1});
            times{m}=temp_str{end-2};
            
            %convert time to minutes and normalize dates
            temp_date=dates(m)-start_date;
                

            temp_hours=times{m}(1:2);
            temp_mins=times{m}(3:4);
            temp_secs=times{m}(5:6);
            
            
            time_mins(m)=str2num(temp_hours)*60+str2num(temp_mins)+...
                str2num(temp_secs)/60+24*60*temp_date-start_time;
            
            temp_str2=strsplit(temp_str{end},'.');
            
            scanner(m)=str2num(temp_str2{1}(1));
            plates{m}=temp_str2{1}(2:end);
            
            %get sga data
            sga_mat{m}=readtable([dir_base exp_name '/' temp_name]);
            
            mat_to_process(:,m)=sga_mat{m}.size;
            
            m=m+1;

        end
        
    end
    
end




plates_to_output={'a','b','e'};


for i=1:length(plates_to_output)
    
    idx_to_use=ismember(plates,plates_to_output{i});
    
    mat_to_output=mat_to_process(:,idx_to_use);
    v_time=time_mins(idx_to_use);

    to_output=table(v_time',mat_to_output');
    writetable(to_output,[dependency_directory 'gitter-data/'...
        condition_names{ismember(plate_names,plates_to_output{i})} '_evo.csv'])
    
end


