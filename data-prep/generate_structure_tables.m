
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/plot'])

tic



load([dependency_directory 'asa_data.mat'])
load([dependency_directory 'neighbor_data.mat'])
load([dependency_directory 'residue_data.mat'])
load([dependency_directory 'secondary_structure_data.mat'])
load([dependency_directory 'superfam_domain_boundaries.mat'])

%make compatible with ismember
residue_mat(cellfun(@isempty,residue_mat))={'NA'};
secondary_mat(cellfun(@isempty,secondary_mat))={'NA'};

load([dependency_directory 'simulated_mutation_positions.mat'])
load([dependency_directory 'simulated_mutation_tstv.mat'])

load([dependency_directory '1K_mutation_positions.mat'])
load([dependency_directory '1K_mutation_af.mat'])



%prepare data matrices for simulated and 1K
asa_mat_1K=nan(size(af_mat_1K));
neighbor_mat_1K=nan(size(af_mat_1K));
residue_mat_1K=nan(size(af_mat_1K));
structure_mat_1K=nan(size(af_mat_1K));
domain_mat_1K=nan(size(af_mat_1K));
ts_mat_1K=nan(size(af_mat_1K));



asa_mat_sim=nan(size(mutation_mat_simulated));
neighbor_mat_sim=nan(size(mutation_mat_simulated));
residue_mat_sim=nan(size(mutation_mat_simulated));
structure_mat_sim=nan(size(mutation_mat_simulated));
domain_mat_sim=nan(size(mutation_mat_simulated));
ts_mat_sim=nan(size(mutation_mat_simulated));


[n_proteins,max_residues]=size(asa_mat);
[~,n_cols]=size(domain_mat);
domain_mat(:,(n_cols+1):max_residues)=0;

[n_rows,~]=size(domain_mat);
domain_mat((n_rows+1):n_proteins,:)=0;

%convert residues to numeric values for storage
for i=1:length(genes_to_use)
    
    i
    
    cols_to_use=~isnan(asa_mat(i,:));
    protein_length(i)=sum(cols_to_use);
    
    temp_asa=asa_mat(i,cols_to_use);
    temp_neighbor=neighbor_mat(i,cols_to_use);
    [~,temp_residue]=residue_types(residue_mat(i,cols_to_use));
    [~,temp_secondary]=structure_types(secondary_mat(i,cols_to_use));
    
    
    pos_1K=mutation_mat_1K(i,:);
    pos_1K=pos_1K(~isnan(pos_1K));
    pos_1K(pos_1K>protein_length(i))=[];

    pos_sim=mutation_mat_simulated(i,:);
    pos_sim=pos_sim(~isnan(pos_sim));
    
    %also build matrices to compare to AF
    asa_mat_1K(i,1:length(pos_1K))=temp_asa(pos_1K);
    asa_mat_sim(i,1:length(pos_sim))=temp_asa(pos_sim);
    
    neighbor_mat_1K(i,1:length(pos_1K))=temp_neighbor(pos_1K);
    neighbor_mat_sim(i,1:length(pos_sim))=temp_neighbor(pos_sim);
    
    residue_mat_1K(i,1:length(pos_1K))=temp_residue(pos_1K);
    structure_mat_1K(i,1:length(pos_1K))=temp_secondary(pos_1K);
        
    residue_mat_sim(i,1:length(pos_sim))=temp_residue(pos_sim);
    structure_mat_sim(i,1:length(pos_sim))=temp_secondary(pos_sim);
    
    
    v_domain_1K=logical(domain_mat(i,pos_1K));
    domain_mat_1K(i,1:length(v_domain_1K))=v_domain_1K;
    
    v_domain_sim=logical(domain_mat(i,pos_sim));
    domain_mat_sim(i,1:length(v_domain_sim))=v_domain_sim;
    
    
    v_ts_1K=logical(ts_mat_simulated(i,pos_1K));
    ts_mat_1K(i,1:length(v_ts_1K))=v_ts_1K;
    
    v_ts_sim=logical(ts_mat_simulated(i,pos_sim));
    ts_mat_sim(i,1:length(v_ts_sim))=v_ts_sim;
    
end





save([dependency_directory 'asa_mat_1K.mat'],'asa_mat_1K')
save([dependency_directory 'neighbor_mat_1K.mat'],'neighbor_mat_1K')
save([dependency_directory 'residue_mat_1K.mat'],'residue_mat_1K')
save([dependency_directory 'structure_mat_1K.mat'],'structure_mat_1K')
save([dependency_directory 'domain_mat_1K.mat'],'domain_mat_1K')
save([dependency_directory 'ts_mat_1K.mat'],'ts_mat_1K')

save([dependency_directory 'asa_mat_sim.mat'],'asa_mat_sim')
save([dependency_directory 'neighbor_mat_sim.mat'],'neighbor_mat_sim')
save([dependency_directory 'residue_mat_sim.mat'],'residue_mat_sim')
save([dependency_directory 'structure_mat_sim.mat'],'structure_mat_sim')
save([dependency_directory 'domain_mat_sim.mat'],'domain_mat_sim')
save([dependency_directory 'ts_mat_sim.mat'],'ts_mat_sim')

save([dependency_directory 'protein_length.mat'],'protein_length')
save([dependency_directory 'gene_names.mat'],'genes_to_use')


toc


