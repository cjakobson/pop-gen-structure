function [] = plot_fpr1_properties(dependency_directory)

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
plot(v_asa_all,'Color',blue)
xlim([55 85])
%highlight mutated residues
v_hits=missense_input.Var1;
scatter(v_hits,(v_asa_all(v_hits)),50,orange,'filled')
ylabel('ASA')
xlabel('residue in Fpr1')
plot(xlim,[30 30],'--k')
axis square


end