function [] = plot_Hs_asa_rare_common(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


hs_data=readtable([dependency_directory 'human-data/misGnomadScrape.csv']);


hs_af=hs_data.Var9/max(hs_data.Var9);

hs_asa=hs_data.sasa;
hs_neighbor=hs_data.neighbors;

[~,hs_structure]=structure_types(hs_data.secondary);


af_thresh=0.001;

temp_idx=hs_af<=af_thresh;

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

hold on
v1=hs_asa(temp_idx);
v1=v1(~isnan(v1));
v2=hs_asa(~temp_idx);
v2=v2(~isnan(v2));

histogram(v1,0:10:300,'Normalization','probability')
histogram(v2,0:10:300,'Normalization','probability')
legend({'rare','common'})
axis square
xlabel('ASA (Ang.^2)')
ylabel('norm. freq')
[h p]=kstest2(v1,v2);
text(250,0.02,num2str(p))




end


