%distance to other alpha carbons from pdb file

clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])

temp_list=dir([dependency_directory 'alphafold-predictions']);


thresh=10;

for i=1:length(temp_list)
    
    if mod(i,100)==0
        i
    end
    
    temp_str=temp_list(i).name;
    
    if length(temp_str)>3
    
        if strcmp(temp_str((end-2):end),'pdb')

            tic
            
            temp_str2=strsplit(temp_str,'-');

            pdb_id=temp_str2{2}
            
            if ~exist([dependency_directory 'neighbor-output/' pdb_id '_' num2str(thresh) 'A.txt'])

                filebase=[dependency_directory 'alphafold-predictions/AF-'];

                file_to_get=[filebase pdb_id '-F1-model_v1.pdb'];

                temp_pdb=readtable(file_to_get,'FileType','text','Delimiter','','ReadVariableNames',0);
                
                clear alpha_mat
                for j=1:height(temp_pdb)

                    temp_str=temp_pdb.Var1{j};

                    if length(temp_str)>3

                        if strcmp(temp_str(1:4),'ATOM')

                            temp_type=temp_str(14:15);

                            if strcmp(temp_type,'CA')

                                temp_res=str2num(temp_str(23:26));

                                temp_coord_x=str2num(temp_str(31:38));
                                temp_coord_y=str2num(temp_str(39:46));
                                temp_coord_z=str2num(temp_str(47:54));

                                alpha_mat(temp_res,:)=[temp_coord_x temp_coord_y temp_coord_z];

                            end

                        end

                    end

                end

                [n_res,~]=size(alpha_mat);

                dist_mat=zeros(n_res);
                for j=1:n_res

                    for k=1:n_res

                        dist_mat(j,k)=sqrt(sum((alpha_mat(j,:)-alpha_mat(k,:)).^2));

                    end

                end

                dist_mat(dist_mat==0)=nan;

                v1=sum(dist_mat<thresh);

                to_ouput=table(v1','VariableNames',{'neighbors'});
                writetable(to_ouput,[dependency_directory 'neighbor-output/' pdb_id '_' num2str(thresh) 'A.txt'])

                toc
                
            end

        end
    
    end
 
end








