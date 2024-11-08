function [] = plot_asa_rare_common_af(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

load([dependency_directory 'residue_mat_1K.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

%test range of af_thresh; make full aa/structure mat for each
%af_thresh=0.01;
af_thresh=10.^[-3:0.1:-1];

for i=1:length(af_thresh)

    temp_idx1=af_mat_1K<=af_thresh(i);

    for j=1:length(aa_labels)

        temp_idx2=residue_mat_1K==j;
        
        for k=1:length(structure_labels)
            
            temp_idx3=structure_mat_1K==k;
            
            v1=reshape(asa_mat_1K(logical(temp_idx1.*temp_idx2.*temp_idx3)),1,[]);
            v2=reshape(asa_mat_1K(logical(~temp_idx1.*temp_idx2.*temp_idx3)),1,[]);
            
            ratio_mat{i}(j,k)=mean(v1,'omitnan')/mean(v2,'omitnan');
            
        end
        
    end

    ratio_mat{i}(isinf(ratio_mat{i}))=nan;
    
    v_ratio=reshape(ratio_mat{i},1,[]);
    v_sd(i)=std(v_ratio,'omitnan')/sqrt(sum(~isnan(v_ratio)));
    v_median(i)=median(v_ratio,'omitnan');

end


%plot(v_median)
hold on
errorbar(v_median,v_sd,'.k')
axis square
title('ASA (Ang.^2)')
ylabel('rare/common')
xlim([0.5 length(af_thresh)+0.5])
ylim([0.85 1.05])
xticks(1:length(af_thresh))
xtickangle(45)
xticklabels(af_thresh)
plot(xlim,[1 1],':r')


end


