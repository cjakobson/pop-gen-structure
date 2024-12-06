function [] = plot_asa_simulation(dependency_directory)

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
%account for double hits in the same pool
missense_input=[missense_input;missense_input(missense_input.Var5==2,:)];

v_asa=dssp_table.sasa(missense_input.Var1);


hold on
histogram(v_asa,0:25:250,'Normalization','probability')
%xlim([0 5])
ylabel('frequency')
xlabel('ASA (Ang.^2)')

n_sims=1e3;
n_clones=length(v_asa);

rng(0)

v_asa_sim=[];
for i=1:n_sims
    
    for j=1:n_clones
        
        rand_idx=ceil(rand*height(mutation_table));
        v_pos(i,j)=mutation_table.residue_number(rand_idx);
        
    end
    
    v_asa_sim=[v_asa_sim; dssp_table.sasa(v_pos(i,:))];
    
end

histogram(v_asa_sim,0:25:250,'Normalization','probability')
axis square
[h p]=kstest2(v_asa,v_asa_sim);
text(100,0.3,num2str(p))
ylim([0 1])
xlabel('ASA (Ang.^2)')




end


