function [] = plot_phi_psi_sim_1K_residue(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'phi_mat_sim.mat'])
load([dependency_directory 'phi_mat_1K.mat'])

load([dependency_directory 'psi_mat_sim.mat'])
load([dependency_directory 'psi_mat_1K.mat'])

load([dependency_directory 'residue_mat_sim.mat'])
load([dependency_directory 'residue_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};


 
figure('units','normalized','outerposition',[0 0 1 1])
p_thresh=1e-4;
%split both by aa and local structure
m=1;
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    for j=1:length(structure_labels)
        
        temp_idx3=structure_mat_sim==j;
        temp_idx4=structure_mat_1K==j;
        
        to_plot{1}=reshape(phi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);
        to_plot{1}(isnan(to_plot{1}))=[];
        %to_plot{3}=reshape(psi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);

        to_plot{2}=reshape(phi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);
        to_plot{2}(isnan(to_plot{2}))=[];
        %to_plot{4}=reshape(psi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);

        [~,p_val1]=kstest2(to_plot{1},to_plot{2});
        %p_val2=ranksum(to_plot{3},to_plot{4});
        
        %only plot if sig diff
        
        if (p_val1<p_thresh)%||(p_val2<p_thresh)
            
            subplot(3,10,m)
            hold on
            histogram(to_plot{1},-180:5:180,'Normalization','Probability')
            histogram(to_plot{2},-180:5:180,'Normalization','Probability')
    %         scatter(v1,v2,1,'k','filled')
    %         scatter(v3,v4,1,'r','filled')
            %easy_box(to_plot)
            axis square
            %xlim([-180 180])
            legend({'1K','sim'})
            xlabel('\phi')
            title([aa_labels{i} ' ' structure_labels{j}])
            
            m=m+1;
        
        end
        
    end

    
end



figure('units','normalized','outerposition',[0 0 1 1])
p_thresh=1e-4;
%split both by aa and local structure
m=1;
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    for j=1:length(structure_labels)
        
        temp_idx3=structure_mat_sim==j;
        temp_idx4=structure_mat_1K==j;
        
        to_plot{1}=reshape(psi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);
        to_plot{1}(isnan(to_plot{1}))=[];
        %to_plot{3}=reshape(psi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);

        to_plot{2}=reshape(psi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);
        to_plot{2}(isnan(to_plot{2}))=[];
        %to_plot{4}=reshape(psi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);

        [~,p_val1]=kstest2(to_plot{1},to_plot{2});
        %p_val2=ranksum(to_plot{3},to_plot{4});
        
        %only plot if sig diff
        
        if (p_val1<p_thresh)%||(p_val2<p_thresh)
            
            subplot(3,10,m)
            hold on
            histogram(to_plot{1},-180:5:180,'Normalization','Probability')
            histogram(to_plot{2},-180:5:180,'Normalization','Probability')
    %         scatter(v1,v2,1,'k','filled')
    %         scatter(v3,v4,1,'r','filled')
            %easy_box(to_plot)
            axis square
            %xlim([-180 180])
            legend({'1K','sim'})
            xlabel('\psi')
            title([aa_labels{i} ' ' structure_labels{j}])
            
            m=m+1;
        
        end
        
    end

    
end





end


