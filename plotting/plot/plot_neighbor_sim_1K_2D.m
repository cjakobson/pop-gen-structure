function [] = plot_neighbor_sim_1K_2D(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'neighbor_mat_sim.mat'])
load([dependency_directory 'neighbor_mat_1K.mat'])

load([dependency_directory 'residue_mat_sim.mat'])
load([dependency_directory 'residue_mat_1K.mat'])

load([dependency_directory 'structure_mat_sim.mat'])
load([dependency_directory 'structure_mat_1K.mat'])

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 
for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    
    temp_neighbor_mat_sim=neighbor_mat_sim(temp_idx1);
    temp_structure_mat_sim=structure_mat_sim(temp_idx1);
    
    temp_neighbor_mat_1K=neighbor_mat_1K(temp_idx2);
    temp_structure_mat_1K=structure_mat_1K(temp_idx2);
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_neighbor_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_neighbor_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
    end
    
end


hold on
imagesc(to_plot,[0.7 1])
yticks(1:length(aa_labels))
yticklabels(aa_labels)
xticks(1:length(structure_labels))
xticklabels(structure_labels)
xtickangle(45)
title('C_\alpha within 10 Ang.')
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


end


