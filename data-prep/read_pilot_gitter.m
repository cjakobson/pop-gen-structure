
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])



exp_name='pilot';

plate_names={'a','b','c','d'};

n_scanners=1;

condition_names={'YPD+rap','NA','NA','NA'};


%get filenames
dir_base=[dependency_directory 'gitter-data/'];

to_read=dir([dir_base exp_name]);


%starting date and time of the exp
start_date=20221023;
start_date2=20221101;
date_offset=8;

start_time=10*60+1+59/60;   %in mins


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
            if dates(m)<start_date2

                temp_date=dates(m)-start_date;
                
            elseif dates(m)>=start_date2
                
                temp_date=dates(m)-start_date2+1+date_offset;
                
            end

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





%just need first plate
idx_to_use=ismember(plates,'a');

mat_to_output=mat_to_process(:,idx_to_use);
v_time=time_mins(idx_to_use);

%split by strain and save
strains={'BY4741','BY4741Dmsh6','BY4743','BY4743DDmsh6'};

%reorganize to 384
a1idxBase=1:2:48;
a2idxBase=2:2:48;
b1idxBase=49:2:96;
b2idxBase=50:2:96;

a1idx=[];
a2idx=[];
b1idx=[];
b2idx=[];

for k=1:16

    a1idx=[a1idx 96*(k-1)+a1idxBase];
    a2idx=[a2idx 96*(k-1)+a2idxBase];
    b1idx=[b1idx 96*(k-1)+b1idxBase];
    b2idx=[b2idx 96*(k-1)+b2idxBase];

end

reorder_mat=nan(size(mat_to_output));

reorder_mat((1:384),:)=mat_to_output(a1idx,:);
reorder_mat(((384+1):(2*384)),:)=mat_to_output(a2idx,:);
reorder_mat(((2*384+1):(3*384)),:)=mat_to_output(b1idx,:);
reorder_mat(((3*384+1):(4*384)),:)=mat_to_output(b2idx,:);


for i=1:length(strains)
    
    idx_to_use=((384*(i-1)+1):(i*384));
    
    to_output=table(v_time',reorder_mat(idx_to_use,:)');
    writetable(to_output,[dependency_directory 'gitter-data/' strains{i} '_pilot.csv'])
    sum(reorder_mat(idx_to_use,end)>1000)
    
end




