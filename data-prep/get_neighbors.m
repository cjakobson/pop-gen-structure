function output_table=get_neighbors(dependency_directory,systematic_name,pdb_id)

file_to_get=[dependency_directory 'neighbor-output/' pdb_id '_10A.txt'];
    
if exist(file_to_get)

    %get neighbors within 10Ang
    temp_neighbors=readtable(file_to_get);
    
    output_table=temp_neighbors;

    save([dependency_directory 'mat-files/' systematic_name '_neighbor_table.mat'],'output_table')
    
end



end