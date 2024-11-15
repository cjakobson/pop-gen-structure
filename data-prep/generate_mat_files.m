
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


input_genome=fastaread([dependency_directory 'orf_coding.fasta']);

systematic_name_dna=cell(length(input_genome),1);
for i=1:height(input_genome)
    
    temp_str=strsplit(input_genome(i).Header,' ');
    systematic_name_dna{i}=temp_str{1};
    
end


for i=1:length(systematic_name_dna)
    
    temp_idx=ismember(systematic_name_pdb,systematic_name_dna{i});
    
    if sum(temp_idx)>0
        temp_pdb_id=pdb_id{temp_idx};

        get_dssp(dependency_directory,systematic_name_dna{i},temp_pdb_id,input_genome(i).Sequence);
        get_neighbors(dependency_directory,systematic_name_dna{i},temp_pdb_id);

        if mod(i,100)==0
            i
        end
    end
    
end



