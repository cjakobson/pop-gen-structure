function [] = plot_mave_cbs_roc(dependency_directory,plot_offset)

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

%mean fitness
unique_residues=unique(v_residue);
for i=1:length(unique_residues)
    
    residue_to_plot(i)=unique_residues(i);
    fitness_to_plot(i)=mean(v_fitness(v_residue==unique_residues(i)));
    
    neighbors_to_plot(i)=neighbor_table.neighbors(unique_residues(i));
    asa_to_plot(i)=dssp_table.sasa(unique_residues(i));
    
end

for i=1:length(v_fitness)

    v_neighbors(i)=neighbor_table.neighbors(v_residue(i));
    v_asa(i)=dssp_table.sasa(v_residue(i));

end


v_asa=v_asa';%asa_to_plot;
v_neighbors=v_neighbors';%neighbors_to_plot;
v_fitness=v_fitness';%fitness_to_plot;





%look up alphamissense predictions
am_data=readtable([dependency_directory 'CBS_alphamissense_AF-P35520-F1-aa-substitutions.csv']);

v_am=nan(length(v_fitness),1);
for i=1:length(v_fitness)

    temp_query=[aminolookup(v_ref{i}) num2str(v_residue(i))...
        aminolookup(v_alt{i})];

    temp_idx=find(ismember(am_data.protein_variant,temp_query));

    if ~isempty(temp_idx)
    
        v_am(i)=am_data.am_pathogenicity(temp_idx);

    end

end

%look up EVE prediction
eve_data=readtable([dependency_directory 'eve_CBS_HUMAN.csv']);

v_eve=nan(length(v_fitness),1);
for i=1:length(v_fitness)

    temp_residue_query=v_residue(i);
    temp_alt_query=aminolookup(v_alt{i});

    temp_idx=find(logical((eve_data.position==temp_residue_query).*...
        ismember(eve_data.mt_aa,temp_alt_query)));

    if ~isempty(temp_idx)
    
        v_eve(i)=eve_data.EVE_scores_ASM(temp_idx);

    end

end






v_bins_asa=[0:20:200 Inf];
v_bins_neighbors=0:3:30;
v_bins_am=0:0.01:1;
v_bins_eve=0:0.01:1;

%threshold for "unfit"

%v_fit_thresh=0.2:0.1:0.4;
%v_unfit_thresh=0.4;
%v_fit_thresh=0.4;
v_unfit_thresh=median(v_fitness,'omitnan');
v_fit_thresh=median(v_fitness,'omitnan');

actual_unfit_idx=v_fitness<=v_unfit_thresh;
actual_fit_idx=v_fitness>=v_fit_thresh;
f_unfit=sum(actual_unfit_idx)/sum(actual_fit_idx);

m=1;
for i=1:length(v_bins_asa)

    for j=1:length(v_bins_neighbors)

        predicted_unfit_idx=logical((v_asa<=v_bins_asa(i)).*...
            (v_neighbors>=v_bins_neighbors(j)));

        tpr(m)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
            sum(actual_unfit_idx);
        fpr(m)=sum(predicted_unfit_idx.*actual_fit_idx)/...
            sum(actual_fit_idx);

        %also PRC
        precision(m)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
                sum(predicted_unfit_idx);
        recall(m)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
                sum(actual_unfit_idx);

        m=m+1;

    end

end

for i=1:length(v_bins_am)

    predicted_unfit_idx=v_am>=v_bins_am(i);

    tpr_am(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
        sum(actual_unfit_idx);
    fpr_am(i)=sum(predicted_unfit_idx.*actual_fit_idx)/...
        sum(actual_fit_idx);

    precision_am(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
            sum(predicted_unfit_idx);
    recall_am(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
            sum(actual_unfit_idx);


    predicted_unfit_idx=v_eve>=v_bins_eve(i);

    tpr_eve(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
        sum(actual_unfit_idx);
    fpr_eve(i)=sum(predicted_unfit_idx.*actual_fit_idx)/...
        sum(actual_fit_idx);

    precision_eve(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
            sum(predicted_unfit_idx);
    recall_eve(i)=sum(predicted_unfit_idx.*actual_unfit_idx)/...
            sum(actual_unfit_idx);

end

subplot(2,4,plot_offset+1)
hold on
scatter(fpr,tpr,10,'k','filled')
scatter(fpr_am,tpr_am,10,'r','filled')
scatter(fpr_eve,tpr_eve,10,'b','filled')

[~,sort_idx]=sort(fpr);
v1=fpr(sort_idx);
v2=tpr(sort_idx);
plot(v1,v2,'-k')
text(0.25,v2(sum(v1<0.25)),num2str(v2(sum(v1<0.25))))

[~,sort_idx]=sort(fpr_am);
v1=fpr_am(sort_idx);
v2=tpr_am(sort_idx);
plot(v1,v2,'-r')
text(0.25,v2(sum(v1<0.25)),num2str(v2(sum(v1<0.25))))

[~,sort_idx]=sort(fpr_eve);
v1=fpr_eve(sort_idx);
v2=tpr_eve(sort_idx);
plot(v1,v2,'-b')
text(0.25,v2(sum(v1<0.25)),num2str(v2(sum(v1<0.25))))


axis square
xlim([0 1])
ylim(xlim)
xlabel('FPR')
ylabel('TPR')
title('human CBS')

plot(xlim,ylim,':r')
text(0.6,0.4,num2str(f_unfit))

legend({'ASA/C_{\alpha}','AlphaMissense','EVE'},'Location','southeast')
%xlabel('C_\alpha within 10 Ang.')


subplot(2,4,plot_offset+2)
hold on
scatter(recall,precision,10,'k','filled')
scatter(recall_am,precision_am,10,'r','filled')
scatter(recall_eve,precision_eve,10,'b','filled')

[~,sort_idx]=sort(recall);
v1=recall(sort_idx);
v2=precision(sort_idx);
plot(v1,v2,'-k')
text(0.75,v2(sum(v1<0.75)),num2str(v2(sum(v1<0.75))))

[~,sort_idx]=sort(recall_am);
v1=recall_am(sort_idx);
v2=precision_am(sort_idx);
plot(v1,v2,'-r')
text(0.75,v2(sum(v1<0.75)),num2str(v2(sum(v1<0.75))))

[~,sort_idx]=sort(recall_eve);
v1=recall_eve(sort_idx);
v2=precision_eve(sort_idx);
plot(v1,v2,'-b')
text(0.75,v2(sum(v1<0.75)),num2str(v2(sum(v1<0.75))))

axis square
xlim([0 1])
ylim(xlim)
xlabel('recall')
ylabel('precision')
title('human CBS')
legend({'ASA/C_{\alpha}','AlphaMissense','EVE'},'Location','southeast')



%repeat with clinvar pathogenic -- just PRC? highly unbalanced





%scatter and do spearman
% subplot(2,4,plot_offset+3)
% hold on
% 
% to_plot1=v_neighbors;
% to_plot2=v_asa;
% 
% temp_idx=v_fitness<=v_unfit_thresh;
% 
% scatter(to_plot1(~temp_idx),to_plot2(~temp_idx),10,'g','filled',...
%     'MarkerFaceAlpha',0.5)
% scatter(to_plot1(temp_idx),to_plot2(temp_idx),10,'m','filled',...
%     'MarkerFaceAlpha',0.5)
% 
% xlabel('C_{\alpha}')
% ylabel('ASA')
% %ylabel('fitness')
% axis square
% %set(gca,'XScale','log')



subplot(2,4,plot_offset+3)
hold on
%v_asa(v_asa==0)=0.1; %jitter off zero
to_plot1=v_neighbors;%./v_asa;
to_plot2=v_am;

temp_idx=v_fitness<=v_unfit_thresh;

scatter(to_plot1(~temp_idx),to_plot2(~temp_idx),5,'g','filled',...
    'MarkerFaceAlpha',0.5)
scatter(to_plot1(temp_idx),to_plot2(temp_idx),5,'m','filled',...
    'MarkerFaceAlpha',0.5)

xlabel('C_{\alpha}')%/ASA')
ylabel('AlphaMissense prediction')
%ylabel('fitness')
axis square
%set(gca,'XScale','log')
[r p]=corr(v_neighbors./v_asa,v_fitness,'rows','complete','type','Spearman');



subplot(2,4,plot_offset+4)
hold on
%v_asa(v_asa==0)=0.1; %jitter off zero
to_plot1=v_asa;
to_plot2=v_am;

temp_idx=v_fitness<=v_unfit_thresh;

scatter(to_plot1(~temp_idx),to_plot2(~temp_idx),5,'g','filled',...
    'MarkerFaceAlpha',0.5)
scatter(to_plot1(temp_idx),to_plot2(temp_idx),5,'m','filled',...
    'MarkerFaceAlpha',0.5)

xlabel('ASA')
ylabel('AlphaMissense prediction')
%ylabel('fitness')
axis square
%set(gca,'XScale','log')
[r p]=corr(v_neighbors./v_asa,v_fitness,'rows','complete','type','Spearman');




%slide fit_thresh and plot AUROC?

