function [] = plot_mave_cbs(dependency_directory,plot_offset)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


cbs_data=readtable([dependency_directory 'mave-datasets/urn_mavedb_00000005-a-4_scores.csv']);

m=1;
%parse mutants
for i=1:height(cbs_data)
    
    temp_str=cbs_data.hgvs_pro{i};
    
    if ~strcmp(temp_str(end),'=')   %skip synonyms
    
        temp_ref=temp_str(3:5);
        temp_alt=temp_str((end-2):end);

        temp_pos=str2num(temp_str(6:(end-3)));

        temp_fitness=cbs_data.score(i);

        v_ref{m}=temp_ref;
        v_alt{m}=temp_alt;
        v_residue(m)=temp_pos;
        v_fitness(m)=temp_fitness;
        m=m+1;
        
    end
    
end


%look up structural parameters
uniprot_id='P35520';

output_table=get_dssp_human(dependency_directory,[],uniprot_id);
dssp_table=output_table;

output_table=get_neighbors_human(dependency_directory,[],uniprot_id);
neighbor_table=output_table;

% %mean fitness
% unique_residues=unique(v_residue);
% for i=1:length(unique_residues)
% 
%     residue_to_plot(i)=unique_residues(i);
%     fitness_to_plot(i)=mean(v_fitness(v_residue==unique_residues(i)));
% 
%     neighbors_to_plot(i)=neighbor_table.neighbors(unique_residues(i));
%     asa_to_plot(i)=dssp_table.sasa(unique_residues(i));
% 
% end

for i=1:length(v_fitness)

    v_neighbors(i)=neighbor_table.neighbors(v_residue(i));
    v_asa(i)=dssp_table.sasa(v_residue(i));

end


%output table to color residues of structures
for i=1:length(dssp_table.sasa)

    v_buffer{i}='';
    v_residue_labels_chimera{i}=[':' num2str(i)];

end

to_output=table(v_buffer',v_residue_labels_chimera',dssp_table.sasa);
writetable(to_output,[dependency_directory 'cbs_sasa.txt'],'Delimiter','\t')

to_output=table(v_buffer',v_residue_labels_chimera',neighbor_table.neighbors);
writetable(to_output,[dependency_directory 'cbs_neighbors.txt'],'Delimiter','\t')


% v1=asa_to_plot;
% v2=fitness_to_plot;

v1=v_asa;
v2=v_fitness;


v_bins=0:20:200;
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
ylim([-0.1 1])
title('human CBS')
for i=1:length(to_plot)
    text(i,0,num2str(length(to_plot{i})))
end



% v1=neighbors_to_plot;
% v2=fitness_to_plot;

v1=v_neighbors;
v2=v_fitness;


v_bins=3:3:30;
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
ylim([-0.1 1])
title('human CBS')
for i=1:length(to_plot)
    text(i,0,num2str(length(to_plot{i})))
end





end


