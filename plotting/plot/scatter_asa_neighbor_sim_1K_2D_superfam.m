function [] = scatter_asa_neighbor_sim_1K_2D_superfam(plot_offset,dependency_directory)


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

to_plot_asa_in_domain=to_plot;

clear to_plot


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
        
    end
    
end

to_plot_asa_outside_domain=to_plot;

clear to_plot






for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    temp_idx5=domain_mat_sim==1;
    temp_idx6=domain_mat_1K==1;
    
    
    temp_neighbor_mat_sim=neighbor_mat_sim(logical(temp_idx1.*temp_idx5));
    temp_structure_mat_sim=structure_mat_sim(logical(temp_idx1.*temp_idx5));
    
    temp_neighbor_mat_1K=neighbor_mat_1K(logical(temp_idx2.*temp_idx6));
    temp_structure_mat_1K=structure_mat_1K(logical(temp_idx2.*temp_idx6));
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_neighbor_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_neighbor_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
    end
    
end

to_plot_neighbor_in_domain=to_plot;



for i=1:length(aa_labels)
    
    temp_idx1=residue_mat_sim==i;
    temp_idx2=residue_mat_1K==i;
    
    temp_idx5=domain_mat_sim==0;
    temp_idx6=domain_mat_1K==0;
    
    
    temp_neighbor_mat_sim=neighbor_mat_sim(logical(temp_idx1.*temp_idx5));
    temp_structure_mat_sim=structure_mat_sim(logical(temp_idx1.*temp_idx5));
    
    temp_neighbor_mat_1K=neighbor_mat_1K(logical(temp_idx2.*temp_idx6));
    temp_structure_mat_1K=structure_mat_1K(logical(temp_idx2.*temp_idx6));
    
    
    for j=1:length(structure_labels)
        
        temp_idx3=temp_structure_mat_sim==j;
        temp_idx4=temp_structure_mat_1K==j;

        to_plot(i,j)=mean(mean(temp_neighbor_mat_1K(temp_idx4),'omitnan'),'omitnan')/...
            mean(mean(temp_neighbor_mat_sim(temp_idx3),'omitnan'),'omitnan');
        
        temp_labels{i,j}=[structure_labels{j} ' ' aa_labels{i}];
        
    end
    
end

to_plot_neighbor_outside_domain=to_plot;





temp_labels=reshape(temp_labels,1,[]);


v1=reshape(to_plot_asa_in_domain,1,[]);
v2=reshape(to_plot_neighbor_in_domain,1,[]);

subplot(2,4,plot_offset+1)
hold on
scatter(v1,v2,25,'k','filled')
xlabel('ASA (Ang.^2) 1K/sim')
ylabel('C_\alpha within 10 Ang.')
title('in domain')
axis square
xlim([1 2.1])
ylim([0.7 1])
%label interesting points
for i=1:length(temp_labels)
    if (v1(i)>1.4)||(v2(i)<0.85)
        text(v1(i),v2(i),temp_labels{i})
    end
end





v1=reshape(to_plot_asa_outside_domain,[],1);
v2=reshape(to_plot_neighbor_outside_domain,[],1);

subplot(2,4,plot_offset+2)
hold on
scatter(v1,v2,25,'k','filled')
xlabel('ASA (Ang.^2) 1K/sim')
ylabel('C_\alpha within 10 Ang.')
title('outside domain')
axis square
xlim([1 2.1])
ylim([0.7 1])
%label interesting points
for i=1:length(temp_labels)
    if (v1(i)>1.4)||(v2(i)<0.85)
        text(v1(i),v2(i),temp_labels{i})
    end
end




end

