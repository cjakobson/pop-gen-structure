function [] = plot_stat3_sasa(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


orfToUse='P40763';  %STAT3

%load([dependency_directory 'human-data/' orfToUse '_neighbor_table.mat'])
%neighbor_table=output_table;
neighbor_table=get_neighbors_human(dependency_directory,'NA',orfToUse);

%load([dependency_directory 'mat-files/' orfToUse '_dssp_table.mat'])
%dssp_table=output_table;
dssp_table=get_dssp_human(dependency_directory,'NA',orfToUse,'NA');

% load([dependency_directory 'mutation-tables/' orfToUse '_mutation_table.mat'])
% mutation_table=output_table;
% mutation_table=mutation_table(mutation_table.is_mis==1,:);


missense_input=readtable([dependency_directory '20260210_stat3_gof_alleles.xlsx']);

for i=1:height(missense_input)

    temp_str=missense_input.ProteinChange{i};

    v_mutated_residue(i)=str2num(temp_str(4:(end-3)));

end



%also gnomad (benign?)
gnomad_input=readtable([dependency_directory 'stat3_gnomAD_v4.1.0_ENSG00000168610_2026_02_10_14_09_37.csv']);
gnomad_input=gnomad_input(ismember(gnomad_input.VEPAnnotation,'missense_variant'),:);

%remove path/likely
gnomad_input=gnomad_input(9:end,:);

for i=1:height(gnomad_input)

    temp_str=gnomad_input.ProteinConsequence{i};

    v_gnomad_residue(i)=str2num(temp_str(6:(end-3)));

end


temp_index=1:height(dssp_table);
temp_index(ismember(temp_index,v_mutated_residue))=[];
to_plot{2}=dssp_table.sasa(temp_index);

to_plot{1}=dssp_table.sasa(v_mutated_residue);

to_plot{3}=dssp_table.sasa(v_gnomad_residue);


hold on
easy_box(to_plot)
xlim([0.5 length(to_plot)+0.5])
ylabel('ASA')
temp_labels={'GoF','all other residues','gnomAD'};
xticks(1:length(temp_labels))
xticklabels(temp_labels)
[h p]=ttest2(to_plot{1},to_plot{2});
text(1.5,200,num2str(p))

end