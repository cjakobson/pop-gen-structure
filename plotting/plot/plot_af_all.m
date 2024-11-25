function [] = plot_af_all(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



load([dependency_directory '1K_mutation_af.mat'])
af_mat_1K(af_mat_1K>0.5)=1-af_mat_1K(af_mat_1K>0.5);

v1=log10(af_mat_1K);
v1=v1(~isnan(v1));
v1=reshape(v1,[],1);



at_data=readtable([dependency_directory 'arabidopsis-data/1001misAnnotated.csv']);

%convert allele counts to MAF
at_af=at_data.AC;
at_af=at_af./max(at_af);
at_af(at_af>0.5)=1-at_af(at_af>0.5);

v2=log10(at_af);
v2=v2(~isnan(v2));



hs_data=readtable([dependency_directory 'human-data/misGnomadScrape.csv']);

hs_af=hs_data.Var9/max(hs_data.Var9);

v3=log10(hs_af);
v3=v3(~isnan(v3));

hold on
histogram(v1,-5:0.2:0,'Normalization','Probability')
histogram(v2,-5:0.2:0,'Normalization','Probability')
histogram(v3,-5:0.2:0,'Normalization','Probability')
set(gca,'YScale','log')
axis square
xlabel('MAF')
ylabel('norm. freq.')
legend({'S.c','A.t.','H.s.'})


end


