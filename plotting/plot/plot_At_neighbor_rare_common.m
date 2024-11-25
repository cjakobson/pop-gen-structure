function [] = plot_At_neighbor_rare_common(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


at_data=readtable([dependency_directory 'arabidopsis-data/1001misAnnotated.csv']);

%convert allele counts to MAF
at_af=at_data.AC;
at_af=at_af./max(at_af);
at_af(at_af>0.5)=1-at_af(at_af>0.5);

at_asa=at_data.sasa;
at_neighbor=at_data.neighbors;

af_thresh=0.001;

temp_idx=at_af<=af_thresh;

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};

hold on
v1=at_neighbor(temp_idx);
v1=v1(~isnan(v1));
v2=at_neighbor(~temp_idx);
v2=v2(~isnan(v2));

histogram(v1,0:2:50,'Normalization','probability')
histogram(v2,0:2:50,'Normalization','probability')
legend({'rare','common'})
axis square
xlabel('C_\alpha within 10 Ang.')
ylabel('norm. freq')
[h p]=kstest2(v1,v2);
text(40,0.02,num2str(p))




end


