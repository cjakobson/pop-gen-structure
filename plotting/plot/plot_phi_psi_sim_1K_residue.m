function [] = plot_phi_psi_sim_1K_residue(dependency_directory,aa_name,structure_name,plot_offset)

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


aa_idx=find(ismember(aa_labels,aa_name));
structure_idx=find(ismember(structure_labels,structure_name));

temp_idx1=residue_mat_sim==aa_idx;
temp_idx2=residue_mat_1K==aa_idx;

temp_idx3=structure_mat_sim==structure_idx;
temp_idx4=structure_mat_1K==structure_idx;


to_plot{1}=reshape(phi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);
to_plot{1}(isnan(to_plot{1}))=[];

to_plot{3}=reshape(psi_mat_1K(logical(temp_idx2.*temp_idx4)),1,[]);
to_plot{3}(isnan(to_plot{1}))=[];


to_plot{2}=reshape(phi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);
to_plot{2}(isnan(to_plot{2}))=[];

to_plot{4}=reshape(psi_mat_sim(logical(temp_idx1.*temp_idx3)),1,[]);
to_plot{4}(isnan(to_plot{2}))=[];


[~,p_val1]=kstest2(to_plot{1},to_plot{2});
[~,p_val2]=kstest2(to_plot{3},to_plot{4});


subplot(2,4,plot_offset+1)
hold on
histogram(to_plot{1},-180:5:180,'Normalization','Probability')
histogram(to_plot{2},-180:5:180,'Normalization','Probability')
axis square
%legend({'1K','sim'})
xlabel('\phi')
ylabel('norm. freq.')
xlim([-180 180])
title([aa_name ' ' structure_name])
text(0,0.05,['p = ' num2str(p_val1)])



subplot(2,4,plot_offset+2)
hold on
histogram(to_plot{3},-180:5:180,'Normalization','Probability')
histogram(to_plot{4},-180:5:180,'Normalization','Probability')
axis square
legend({'1K','sim'},'Location','northwest')
xlabel('\psi')
ylabel('norm. freq.')
xlim([-180 180])
title([aa_name ' ' structure_name])
text(-150,0.03,['p = ' num2str(p_val2)])



subplot(2,4,plot_offset+3)
hold on
h1=histogram2(to_plot{2},to_plot{4},-180:10:180,-180:10:180,...
    'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability');
title('sim')
xlabel('\phi')
ylabel('\psi')
xlim([-180 180])
ylim([-180 180])
axis square



subplot(2,4,plot_offset+4)
hold on
h2=histogram2(to_plot{1},to_plot{3},-180:10:180,-180:10:180,...
    'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability');
title('1K')
xlabel('\phi')
ylabel('\psi')
xlim([-180 180])
ylim([-180 180])
axis square








end


