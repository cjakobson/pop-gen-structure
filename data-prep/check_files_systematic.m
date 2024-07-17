
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


temp_list=pdb_id;
%check for presence of DSSP and neighbor files
m=1;
n=1;
o=1;
thresh=10;
for i=1:length(temp_list)
    
    temp_pdb_id=temp_list{i};
            
    if ~exist([dependency_directory 'neighbor-output/' temp_pdb_id '_' num2str(thresh) 'A.txt'])

        missing_neighbors{m}=systematic_name_pdb{i};
        m=m+1;

    end

    if ~exist([dependency_directory 'dssp-output/AF-' temp_pdb_id '-F1-model_v1.pdbdssp.txt'])

        missing_dssp{n}=systematic_name_pdb{i};
        n=n+1;

    end

    if ~exist([dependency_directory 'alphafold-predictions/AF-' temp_pdb_id '-F1-model_v1.pdb'])

        missing_af{o}=systematic_name_pdb{i};
        o=o+1;

    end
            
end

m-1
n-1
o-1

to_output=table(missing_af','VariableNames',{'ORF'});
writetable(to_output,[dependency_directory 'missing_genes.txt'])


