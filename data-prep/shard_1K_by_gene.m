
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])

%load genes to use
load([dependency_directory 'asa_data.mat'])

load([dependency_directory '1K_data_annotated.mat'])

%fill empty cells
gene(cellfun(@isempty,gene))={'NA'};
type(cellfun(@isempty,type))={'NA'};

tic

temp_idx1=ismember(type,'missense_variant');
for i=1:length(genes_to_use)
    
    query_gene=genes_to_use{i};
    
    temp_idx2=ismember(gene,query_gene);
    
    temp_idx=logical(temp_idx1.*temp_idx2);
    
    temp_encoded=proteinEncoded(temp_idx);
    
    v_residue=nan(sum(temp_idx),1);
    v_af=af(temp_idx);
    v_ref=cell(length(v_residue),1);
    v_alt=cell(length(v_residue),1);
    
    for j=1:length(temp_encoded)
        
        v_ref{j}=temp_encoded{j}(3:5);
        v_alt{j}=temp_encoded{j}((end-2):end);
        
        v_residue(j)=str2num(temp_encoded{j}(6:(end-3)));
        
    end
    
    to_output=table(v_residue,v_ref,v_alt,v_af);
    
    save([dependency_directory '1K-data-by-gene/' query_gene '_1K_data.mat'],'to_output')
    
    n_mutations(i)=length(v_residue);
    
end


%make matrix of mutation positions in same order as simulation matrices
max(n_mutations)

mutation_mat_1K=nan(length(genes_to_use),max(n_mutations));
af_mat_1K=nan(length(genes_to_use),max(n_mutations));

for i=1:length(genes_to_use)
    
    query_gene=genes_to_use{i};
    
    load([dependency_directory '1K-data-by-gene/' query_gene '_1K_data.mat'])
    
    if height(to_output)>0
        
        mutation_mat_1K(i,1:height(to_output))=to_output.v_residue;
        %also allele frequencies
        af_mat_1K(i,1:height(to_output))=to_output.v_af;
        
    end
    
end

save([dependency_directory '1K_mutation_positions.mat'],'mutation_mat_1K')
save([dependency_directory '1K_mutation_af.mat'],'af_mat_1K')

toc


tic

%do the same for all possible mutations
for i=1:length(genes_to_use)
    
    query_gene=genes_to_use{i};
    
    load([dependency_directory 'mutation-tables/' query_gene '_mutation_table.mat'])
    
    n_possible(i)=sum(output_table.is_mis);
    
end


max(n_possible)




mutation_mat_simulated=nan(length(genes_to_use),max(n_possible));
ts_mat_simulated=nan(length(genes_to_use),max(n_possible));

for i=1:length(genes_to_use)
    
    query_gene=genes_to_use{i};
    
    load([dependency_directory 'mutation-tables/' query_gene '_mutation_table.mat'])
    
    temp_idx=output_table.is_mis==1;
    
    mutation_mat_simulated(i,1:sum(temp_idx))=output_table.residue_number(temp_idx);
    
    %also Ts/Tv
    ts_mat_simulated(i,1:sum(temp_idx))=output_table.is_ts(temp_idx);
    
    
end


save([dependency_directory 'simulated_mutation_positions.mat'],'mutation_mat_simulated')
save([dependency_directory 'simulated_mutation_tstv.mat'],'ts_mat_simulated')



%superfam domain boundaries
domain_input=readtable([dependency_directory 'Saccharomyces_cerevisiae_SUPERFAMILY_domains.txt']);

for i=1:length(genes_to_use)
    
    query_gene=genes_to_use{i};
    
    temp_idx=find(ismember(domain_input.SequenceID,query_gene));
    
    v_temp=[];
    
    if length(temp_idx)>0
    
        for j=1:length(temp_idx)

            temp_boundaries=strsplit(domain_input.RegionOfAssignment{temp_idx(j)},',');

            for k=1:length(temp_boundaries)

                temp_str=strsplit(temp_boundaries{k},'-');

                left_bound=str2num(temp_str{1});
                right_bound=str2num(temp_str{2});

                v_temp(left_bound:right_bound)=1;

            end

        end

        domain_mat(i,1:length(v_temp))=v_temp;
        
    end
    
end


save([dependency_directory 'superfam_domain_boundaries.mat'],'domain_mat')

toc


