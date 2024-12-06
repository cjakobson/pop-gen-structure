function [] = scatter_evo_asa_neighbors(dependency_directory)

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



%all possible missense
v_asa_all=dssp_table.sasa;
v_neighbors_all=neighbor_table.neighbors;

v_asa=dssp_table.sasa(missense_input.Var1);
v_neighbors=neighbor_table.neighbors(missense_input.Var1);

hold on
scatter(v_asa_all,v_neighbors_all,25,'k','filled')
scatter(v_asa,v_neighbors,50,'r','filled')
axis square
xlim([-10 250])
ylim([0 30])
ylabel('C_\alpha within 10 Ang.')
xlabel('ASA (Ang.^2)')

asa_thresh=30;
neighbor_thresh=15;

n_possible=sum((v_asa_all<=asa_thresh).*(v_neighbors_all>=neighbor_thresh));
possible_residues=find((v_asa_all<=asa_thresh).*(v_neighbors_all>=neighbor_thresh));
n_found=sum(ismember(possible_residues,missense_input.Var1));
n_clones=sum(ismember(missense_input.Var1,possible_residues));
%n_found=sum((v_asa<=asa_thresh).*(v_neighbors>=neighbor_thresh));

text(50,30,['n possible = ' num2str(n_possible)])
text(50,28,['n found = ' num2str(n_found)])
text(50,26,['n clones = ' num2str(n_clones)])



end