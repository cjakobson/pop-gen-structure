function [] = plot_multi_hit_simulation(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


orfToUse='YNL135C';

load([dependency_directory 'mat-files/' orfToUse '_neighbor_table.mat'])
neighbor_table=output_table;

load([dependency_directory 'mat-files/' orfToUse '_dssp_table.mat'])
dssp_table=output_table;

load([dependency_directory 'mutation-tables/' orfToUse '_mutation_table.mat'])
mutation_table=output_table;
mutation_table=mutation_table(mutation_table.is_mis==1,:);


missense_input=readtable([dependency_directory 'FPR1_primordium_missense.csv']);
lof_input=readtable([dependency_directory 'FPR1_primordium_lof.csv']);

unique_missense=unique(missense_input.Var2);

m=1;
for i=1:length(unique_missense)
    
    v_hits(m)=sum(missense_input.Var5(ismember(missense_input.Var2,unique_missense{i})));
    m=m+1;
    
end

unique_lof=unique(lof_input.Var2);

for i=1:length(unique_lof)
    
    v_hits(m)=sum(lof_input.Var5(ismember(lof_input.Var2,unique_lof{i})));
    m=m+1;
    
end

hold on
histogram(v_hits,'Normalization','probability')
xlim([0 5])
ylabel('frequency')
xlabel('hits per mutation')

n_sims=1e3;
n_clones=sum(v_hits);

rng(0)

v_hits_sim=[];
for i=1:n_sims
    
    for j=1:n_clones
        
        rand_idx=ceil(rand*height(mutation_table));
        v_pos(i,j)=mutation_table.residue_number(rand_idx);
        
    end
    
    unique_pos=unique(v_pos(i,:));
    for j=1:length(unique_pos)
        
        sim_frequency{i}(j)=sum(v_pos(i,:)==unique_pos(j));
        
    end
    
    v_hits_sim=[v_hits_sim sim_frequency{i}];
    
end

histogram(v_hits_sim,'Normalization','probability')
axis square
ylim([0 1])



end


