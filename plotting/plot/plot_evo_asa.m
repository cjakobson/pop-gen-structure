function [] = plot_evo_asa_neighbors(dependency_directory)

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

to_plot{1}=v_asa;

%all possible missense
to_plot{2}=dssp_table.sasa(mutation_table.residue_number);

hold on
easy_box(to_plot)
[h p]=ttest2(to_plot{1},to_plot{2});
ylim([-10 250])
text(1.5,200,num2str(p))
xticklabels({'evo','all'})
ylabel('ASA (Ang.^2)')

end


