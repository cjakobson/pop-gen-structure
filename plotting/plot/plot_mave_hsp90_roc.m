function [] = plot_mave_hsp90_roc(dependency_directory,plot_offset)

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
v2=neighbors_to_plot;
v3=fitness_to_plot;

v_bins_asa=0:30:180;
v_bins_neighbors=0:5:30;

%threshold for "unfit"

%v_fit_thresh=-1:0.25:-0.25
v_unfit_thresh=-0.2;
v_fit_thresh=-0.1;

for k=1:length(v_fit_thresh)

    actual_unfit_idx=v3<v_unfit_thresh(k);
    actual_fit_idx=v3>v_fit_thresh(k);
    f_unfit=sum(actual_unfit_idx)/sum(actual_fit_idx)
    text(0.6,0.1,num2str(f_unfit))
    
    m=1;
    for i=1:length(v_bins_asa)
    
        for j=1:length(v_bins_neighbors)
    
            predicted_unfit_idx=logical((v1<v_bins_asa(i)).*...
                (v2>v_bins_neighbors(j)));
    
            tpr(m)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
                sum(actual_unfit_idx);
            fpr(m)=sum(predicted_unfit_idx.*actual_fit_idx)/...
                sum(actual_fit_idx);
            m=m+1;
    
        end
    
    end
    
    hold on
    scatter(fpr,tpr,10*k,'k','filled')
    %plot(fpr,tpr,'-k')
    axis square
    xlim([0 1])
    ylim(xlim)
    xlabel('FPR')
    ylabel('TPR')
    title('yeast Hsp90')
    
    plot(xlim,ylim,':r')
    %xlabel('C_\alpha within 10 Ang.')
    
    %xlabel('ASA (Ang.^2)')

end




end


