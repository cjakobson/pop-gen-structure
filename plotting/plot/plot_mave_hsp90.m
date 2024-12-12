function [] = plot_mave_hsp90(dependency_directory,plot_offset)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


hsp90_data=readtable([dependency_directory 'mave-datasets/urn_mavedb_00000074-a-1_scores.csv']);

m=1;
%parse mutants
for i=1:height(hsp90_data)
    
    temp_str=hsp90_data.hgvs_pro{i};
    
    if ~strcmp(temp_str(end),'=')   %skip synonyms
    
        temp_ref=temp_str(3:5);
        temp_alt=temp_str((end-2):end);

        temp_pos=str2num(temp_str(6:(end-3)));

        temp_fitness=hsp90_data.score(i);

        v_ref{m}=temp_ref;
        v_alt{m}=temp_alt;
        v_residue(m)=temp_pos;
        v_fitness(m)=temp_fitness;
        m=m+1;
        
    end
    
end


%look up structural parameters
systematic_name='YPL240C';

load([dependency_directory 'mat-files/' systematic_name '_dssp_table.mat'])
dssp_table=output_table;

load([dependency_directory 'mat-files/' systematic_name '_neighbor_table.mat'])
neighbor_table=output_table;

%mean fitness
unique_residues=unique(v_residue);
for i=1:length(unique_residues)
    
    residue_to_plot(i)=unique_residues(i);
    fitness_to_plot(i)=mean(v_fitness(v_residue==unique_residues(i)));
    
    neighbors_to_plot(i)=neighbor_table.neighbors(unique_residues(i));
    asa_to_plot(i)=dssp_table.sasa(unique_residues(i));
    
end


v1=asa_to_plot;
v2=fitness_to_plot;

v_bins=0:30:180;
clear to_plot
for i=2:length(v_bins)
    temp_idx=logical((v1>=v_bins(i-1)).*(v1<v_bins(i)));
    to_plot{i-1}=v2(temp_idx);
end

subplot(2,4,plot_offset+1)
easy_box(to_plot)
axis square
xlabel('ASA (Ang.^2)')
ylabel('fitness')
xticks(1:length(v_bins(2:end)))
xticklabels(v_bins(2:end))
ylim([-1.5 0.1])
title('yeast Hsp90')



v1=neighbors_to_plot;
v2=fitness_to_plot;

v_bins=0:5:30;
clear to_plot
for i=2:length(v_bins)
    temp_idx=logical((v1>=v_bins(i-1)).*(v1<v_bins(i)));
    to_plot{i-1}=v2(temp_idx);
end

subplot(2,4,plot_offset+2)
easy_box(to_plot)
axis square
xlabel('C_\alpha within 10 Ang.')
ylabel('fitness')
xticks(1:length(v_bins(2:end)))
xticklabels(v_bins(2:end))
ylim([-1.5 0.1])
title('yeast Hsp90')




end

