
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/plot'])

tic

load([dependency_directory 'asa_mat_1K.mat'])

sc_n_muts=sum(~isnan(asa_mat_1K),2);


at_data=readtable([dependency_directory 'arabidopsis-data/1001misAnnotated.csv']);

at_v_gene=at_data.uniprotName;
at_genes=unique(at_v_gene);
length(at_genes)
at_n_muts=nan(length(at_genes),1);
for i=1:length(at_genes)

    if mod(i,1000)==0
        i
    end

    at_n_muts(i)=sum(ismember(at_v_gene,at_genes{i}));

end


hs_data=readtable([dependency_directory 'human-data/misGnomadScrape.csv']);

hs_v_gene=hs_data.Var14;
hs_genes=unique(hs_v_gene);
length(hs_genes)
hs_n_muts=nan(length(hs_genes),1);
for i=1:length(hs_genes)

    if mod(i,1000)==0
        i
    end

    hs_n_muts(i)=sum(ismember(hs_v_gene,hs_genes{i}));

end




save([dependency_directory 'mutations_per_gene.mat'],'sc_n_muts',...
    'at_n_muts','hs_n_muts')



toc





