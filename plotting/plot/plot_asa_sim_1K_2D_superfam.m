function [] = plot_asa_sim_1K_2D_superfam(plot_offset,dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])

load([dependency_directory 'residue_mat_sim.mat'])
load([dependency_directory 'residue_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

load([dependency_directory 'domain_mat_sim.mat'])
load([dependency_directory 'domain_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    temp_idx5=domain_mat_sim==1;
    temp_idx6=domain_mat_1K==1;
    
    temp_asa_mat_sim=asa_mat_sim(logical(temp_idx1.*temp_idx5));
    temp_structure_mat_sim=structure_mat_sim(logical(temp_idx1.*temp_idx5));
    
    temp_asa_mat_1K=asa_mat_1K(logical(temp_idx2.*temp_idx6));
    temp_structure_mat_1K=structure_mat_1K(logical(temp_idx2.*temp_idx6));
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_asa_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_asa_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
    end
    
end


subplot(2,8,plot_offset+1)
hold on
imagesc(to_plot,[1 1.7])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('ASA (Ang.^2) 1K/sim in domain')
xlim([0.5 length(structure_labels)+.5])
ylim([0.5 length(aa_labels)+.5])

m = size(get(gcf,'colormap'),1);
%red to blue colormap
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
%colormap(flipud(c))
colormap(c)
colorbar

to_plot_in_domain=to_plot;


m=1;
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    temp_idx5=domain_mat_sim==0;
    temp_idx6=domain_mat_1K==0;
    
    temp_asa_mat_sim=asa_mat_sim(logical(temp_idx1.*temp_idx5));
    temp_structure_mat_sim=structure_mat_sim(logical(temp_idx1.*temp_idx5));
    
    temp_asa_mat_1K=asa_mat_1K(logical(temp_idx2.*temp_idx6));
    temp_structure_mat_1K=structure_mat_1K(logical(temp_idx2.*temp_idx6));
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_asa_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_asa_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
        temp_labels{m}=[structure_labels{j} ' ' aa_labels{i}];
        m=m+1;
        
    end
    
end


subplot(2,8,plot_offset+2)
hold on
imagesc(to_plot,[1 1.7])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('ASA (Ang.^2) 1K/sim outside domain')
xlim([0.5 length(structure_labels)+.5])
ylim([0.5 length(aa_labels)+.5])
colormap(c)
colorbar


to_plot_outside_domain=to_plot;
% 
% 
% %scatter within vs outside
% v1=reshape(to_plot_outside_domain,1,[]);
% v2=reshape(to_plot_in_domain,1,[]);
% 
% subplot(2,4,plot_offset/2+2)
% scatter(v1,v2,10,'k','filled')
% for i=1:length(temp_labels)
%     if abs(v1(i)-v2(i))>0.2
%         text(v1(i),v2(i),temp_labels{i})
%     end
% end
% axis square
% xlabel('outside domain')
% ylabel('in domain')

end


