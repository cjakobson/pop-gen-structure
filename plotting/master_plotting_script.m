
clear

tic


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)


%filebase='/Users/cjakobson/';
filebase='/Users/christopherjakobson/';

figure_counter=1;

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-output/'];


addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/plot'])




%Figure 1
figure('units','normalized','outerposition',[0 0 1 1])

%A
%scheme


%ASA and neighbors for simulated mutations
%B
subplot(2,4,1)
plot_asa_sim(dependency_directory)

%C
subplot(2,4,2)
plot_neighbor_sim(dependency_directory)



%D
%plot relative to 1K
subplot(2,4,3)
plot_asa_sim_1K(dependency_directory)

%E
subplot(2,4,4)
plot_neighbor_sim_1K(dependency_directory)



%break down by secondary structure and AA
%F
subplot(2,4,5)
plot_sim_1K_structure(dependency_directory)


%G
subplot(2,2,4)
plot_sim_1K_residue(dependency_directory)





set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_1_1'],'-dsvg','-r0')
print([output_directory 'Figure_1_1'],'-djpeg','-r300')



%Figure S1

%raw bar charts by structure and residue
figure('units','normalized','outerposition',[0 0 1 1])

%A
subplot(2,4,1)
plot_asa_sim_raw_structure(dependency_directory)

%B
subplot(2,4,2)
plot_neighbor_sim_raw_structure(dependency_directory)

%C
subplot(2,4,3)
plot_asa_sim_raw_residue(dependency_directory)

%D
subplot(2,4,4)
plot_neighbor_sim_raw_residue(dependency_directory)



%Ts/Tv




set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_S1_1'],'-dsvg','-r0')
print([output_directory 'Figure_S1_1'],'-djpeg','-r300')





%Figure 2
%gene-level analysis
figure('units','normalized','outerposition',[0 0 1 1])


%gene-level strength of selection is correlated across metrics
%A
subplot(2,4,1)
scatter_asa_neighbor_gene(dependency_directory,output_directory)

%B
subplot(2,4,2)
scatter_asa_dN_gene(dependency_directory)

%C
subplot(2,4,3)
scatter_neighbor_dN_gene(dependency_directory)


%GO-term analysis -- highlight metabolic process genes on A? boxplot of
%split by GO term?



%E
%scheme


%split genes into young/ancient
%F
subplot(2,4,4)
plot_asa_sim_1K_structure_age(dependency_directory)

%G
subplot(2,4,5)
plot_neighbor_sim_1K_structure_age(dependency_directory)





set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_2_1'],'-dsvg','-r0')
print([output_directory 'Figure_2_1'],'-djpeg','-r300')





figure('units','normalized','outerposition',[0 0 1 1])

%2-D plots by AA/struct
%subplot(2,8,1)
%plot_asa_sim_1K_2D(dependency_directory)

%subplot(2,8,2)
%plot_neighbor_sim_1K_2D(dependency_directory)




%scatter by combo
subplot(2,4,1)
scatter_asa_neighbor_sim_1K_2D(dependency_directory)



%metrics by within domain/not
subplot(2,8,3)
plot_asa_sim_1K_superfam(dependency_directory)

subplot(2,8,4)
plot_neighbor_sim_1K_superfam(dependency_directory)




%break down phi/psi for unstr G specifically

%put 2D plots in supp
plot_phi_psi_sim_1K_residue(dependency_directory,'G','unstr.',4)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_2_2'],'-dsvg','-r0')
print([output_directory 'Figure_2_2'],'-djpeg','-r300')





%repeat 2D analysis this way
%plot_asa_sim_1K_2D_superfam(8,dependency_directory)

%plot_neighbor_sim_1K_2D_superfam(10,dependency_directory)




%scatter_asa_neighbor_sim_1K_2D_superfam(6,dependency_directory)




%Figure S2
figure('units','normalized','outerposition',[0 0 1 1])



%also by residue
%A
subplot(2,4,1)
plot_asa_sim_1K_residue_age(dependency_directory)

%B
subplot(2,4,2)
plot_neighbor_sim_1K_residue_age(dependency_directory)

%C
subplot(2,4,3)
plot_dN_sim_1K_residue_age(dependency_directory)

%D
subplot(2,4,4)
plot_dN_sim_1K_structure_age(dependency_directory)



%put proline in the supp for reference
plot_phi_psi_sim_1K_residue(dependency_directory,'P','unstr.',4)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_S2_1'],'-dsvg','-r0')
print([output_directory 'Figure_S2_1'],'-djpeg','-r300')





%Figure 3
%rare/common in 1K
figure('units','normalized','outerposition',[0 0 1 1])


%A
%scheme

%B
subplot(2,4,1)
plot_asa_rare_common(dependency_directory)

%C
subplot(2,4,2)
plot_neighbor_rare_common(dependency_directory)

%D
%merge these
subplot(2,4,3)
plot_rare_common_structure(dependency_directory)



%E
subplot(2,2,3)
plot_rare_common_residue(dependency_directory)




%slide af_thresh

%F
subplot(2,4,7)
plot_neighbor_rare_common_af(dependency_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_3_1'],'-dsvg','-r0')
print([output_directory 'Figure_3_1'],'-djpeg','-r300')



%niche enrichments
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,4,1)





%Figure S3
figure('units','normalized','outerposition',[0 0 1 1])

%A
subplot(2,4,1)
plot_asa_rare_common_af(dependency_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'Figure_S3_1'],'-dsvg','-r0')
print([output_directory 'Figure_S3_1'],'-djpeg','-r300')



toc


