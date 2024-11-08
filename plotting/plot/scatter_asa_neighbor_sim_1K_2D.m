function [] = scatter_asa_neighbor_sim_1K_2D(dependency_directory)


blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'asa_mat_sim.mat'])
load([dependency_directory 'asa_mat_1K.mat'])


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
    
    
    temp_asa_mat_sim=asa_mat_sim(temp_idx1);
    temp_structure_mat_sim=structure_mat_sim(temp_idx1);
    
    temp_asa_mat_1K=asa_mat_1K(temp_idx2);
    temp_structure_mat_1K=structure_mat_1K(temp_idx2);
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_asa_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_asa_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
    end
    
end

to_plot_asa=to_plot;

clear to_plot


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

to_plot_neighbor=to_plot;

v1=reshape(to_plot_asa,1,[]);
v2=reshape(to_plot_neighbor,1,[]);

hold on
scatter(v1,v2,25,'k','filled')
xlabel('ASA (Ang.^2) 1K/sim')
ylabel('C_\alpha within 10 Ang.')
axis square
xlim([1 1.7])
ylim([0.7 1])



end


