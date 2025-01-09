function [] = plot_1K_roc(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

load([dependency_directory 'residue_mat_1K.mat'])
load([dependency_directory 'asa_mat_1K.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

for n=1:length(aa_labels)

    temp_idx=residue_mat_1K==n;

    v1=reshape(asa_mat_1K(temp_idx),[],1);
    v2=reshape(neighbor_mat_1K(temp_idx),[],1);
    v3=reshape(af_mat_1K(temp_idx),[],1);

    
    v_bins_asa=0:30:180;
    v_bins_neighbors=0:5:30;
    
    %threshold for "unfit"
    common_thresh=0.1;
    rare_thresh=0.0005;
    
    for k=1:length(common_thresh)
    
        actual_unfit_idx=v3<rare_thresh(k);
        actual_fit_idx=v3>common_thresh(k);
        f_unfit=sum(actual_unfit_idx)/sum(actual_fit_idx);
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
        
        subplot(4,5,n)
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


end


