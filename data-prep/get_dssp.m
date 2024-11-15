function output_table=get_dssp(dependency_directory,systematic_name,pdb_id,dna_sequence)

fileToGet=[dependency_directory 'dssp-output/AF-' pdb_id '-F1-model_v1.pdbdssp.txt'];
    
if exist(fileToGet)

    %get DSSP
    temp_dssp=readtable(fileToGet,'ReadVariableNames',0,'delimiter',' A ');

    secondary=cell(height(temp_dssp),1);
    residue=cell(height(temp_dssp),1);
    no_struct=zeros(height(temp_dssp),1);
    sasa=zeros(height(temp_dssp),1);
    phi=zeros(height(temp_dssp),1);
    psi=zeros(height(temp_dssp),1);
    pos_in_orf=zeros(height(temp_dssp),1);
    
    for j=1:height(temp_dssp)
        
        temp_str=temp_dssp.Var2{j};
        
        if length(temp_str)>=4
            
            residue{j}=temp_str(1);
            secondary{j}=temp_str(4);
            
            if strcmp(temp_str(4),' ')
                
                no_struct(j)=1;
                
            end
            
            if length(temp_str)>=25
                
                sasa(j)=str2num(temp_str(23:25));
                phi(j)=str2num(temp_str(91:96));
                psi(j)=str2num(temp_str(97:102));
                
            end
            
        end
        
        pos_in_orf(j)=j/height(temp_dssp);
        
    end

    %calculate lengths of runs in secondary structure (to
    %find positions in 'metahelix', etc
    temp_is_run=zeros(length(secondary),1);
    for j=2:length(secondary)
        if strcmp(secondary{j},secondary{j-1})
            temp_is_run(j-1)=1;
        end
    end

    temp_pos_in_run=zeros(length(secondary),1);
    temp_counter=0;
    for j=1:length(temp_pos_in_run)
        if temp_is_run(j)>0
            temp_counter=temp_counter+temp_is_run(j);
            temp_pos_in_run(j)=temp_counter;
        else
            temp_counter=0;
        end
    end

    temp_lengthof_run=temp_pos_in_run;
    for j=1:(length(temp_lengthof_run)-1)
        temp_idx=length(temp_lengthof_run)-(j-1);
        if temp_pos_in_run(temp_idx)==0
            temp_lengthof_run(temp_idx)=0;
        elseif temp_lengthof_run(temp_idx)>temp_lengthof_run(temp_idx-1)
            temp_lengthof_run(temp_idx-1)=temp_lengthof_run(temp_idx);
        end
    end
    
    for j=1:length(secondary)
        if strcmp(secondary{j},' ')
            secondary{j}='U';    %use U for unstructured (as opposed to blank which is missing data)
        end
    end

    pos_in_run=temp_pos_in_run;
    length_of_run=temp_lengthof_run;
    
    output_table=table(residue,secondary,sasa,phi,psi,pos_in_orf,pos_in_run,length_of_run);

    save([dependency_directory 'mat-files/' systematic_name '_dssp_table.mat'],'output_table','dna_sequence')
    
end

    





end