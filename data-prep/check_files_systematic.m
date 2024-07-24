
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])


pdb_dictionary=tdfread([dependency_directory 'uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.07.07-20.59.31.22.tsv']);

pdb_id=cell(length(pdb_dictionary.Entry),1);
systematic_name_pdb=pdb_id;
m=1;
for i=1:length(pdb_dictionary.Entry)
    
    temp_str=pdb_dictionary.Entry(i,:);
    temp_str(temp_str==' ')=[];
    pdb_id{m}=temp_str;
    
    temp_str=pdb_dictionary.Gene_Names_0x28ordered_locus0x29(i,:);
    temp_str(temp_str==' ')=[];
    systematic_name_pdb{m}=temp_str;
    m=m+1;
    
end


tic

temp_list=pdb_id;
%check for presence of DSSP and neighbor files, etc
%also check lengths of sequences
seq_length=nan(length(temp_list),5);
thresh=10;
for i=1:length(temp_list)
    
    temp_pdb_id=temp_list{i};
            
    file_to_get{1}=[dependency_directory 'neighbor-output/' temp_pdb_id '_' num2str(thresh) 'A.txt'];
    file_to_get{2}=[dependency_directory 'dssp-output/AF-' temp_pdb_id '-F1-model_v1.pdbdssp.txt'];
    file_to_get{3}=[dependency_directory 'mat-files/' systematic_name_pdb{i} '_neighbor_table.mat'];
    file_to_get{4}=[dependency_directory 'mat-files/' systematic_name_pdb{i} '_dssp_table.mat'];
    file_to_get{5}=[dependency_directory 'mutation-tables/' systematic_name_pdb{i} '_mutation_table.mat'];
    
    if exist(file_to_get{1})
        
        temp_table=readtable(file_to_get{1});
        seq_length(i,1)=height(temp_table);

    end
    
    if exist(file_to_get{2})
        
        temp_table=readtable(file_to_get{2});
        seq_length(i,2)=height(temp_table)-21;  %header is 21 lines

    end
    
    if exist(file_to_get{3})
        
        load(file_to_get{3})
        seq_length(i,3)=height(output_table);

    end
    
    if exist(file_to_get{4})
        
        load(file_to_get{4})
        seq_length(i,4)=height(output_table);

    end
    
    if exist(file_to_get{5})
        
        load(file_to_get{5})
        seq_length(i,5)=height(output_table)/9; %3 nt per aa, 3 muts per nt

    end
    
            
end

toc


%missing data
sum(isnan(seq_length))
% to_output=table(missing_dssp_mat','VariableNames',{'ORF'});
% writetable(to_output,[dependency_directory 'missing_genes.txt'])



%check whether pdb and genome protein lengths, etc agree
sum(seq_length(:,3)~=seq_length(:,4))
sum(seq_length(:,4)~=seq_length(:,5))

%check which genes have no missing data
sum(sum(isnan(seq_length),2)==0)

%check which genes have no length mismatches




%how to format ASA and neighbor data? genes as rows and residues as
%columns? then references this for 1K data


