
clear

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/pop-gen-structure-dependencies/'];

addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/plot'])

tic


ecotype_data=readtable([dependency_directory '1002_genomes_ecotypes.txt'],...
    'NumHeaderLines',0);


load([dependency_directory '1002_data_common.mat'])


strain_names=ecotype_data.StandardID;

%match strain names
for i=1:length(strainString)
    
    if length(strainString{i})>4
        
        strain_query=[strainString{i}((end-2):end) '*'];
        
    else
        
        strain_query=strainString{i};
        
    end

    temp_idx=find(ismember(strain_names,strain_query));

    if ~isempty(temp_idx)
            
        v_niche{i}=ecotype_data.Ecotype{temp_idx};

    else

        v_niche{i}='NA';

    end

end
niche_names=unique(v_niche);

sum(ismember(v_niche,'NA'));



%convert to numerical
for i=1:length(niche_names)

    v_niche_code(ismember(v_niche,niche_names{i}))=i;

end

%calculate coherence per-mutation and plot against MAF
niche_mat_ref=nan(length(chr),length(niche_names));
niche_mat_alt=nan(length(chr),length(niche_names));
for i=1:length(chr)

    v_temp=minGenotype(i,:);

    ref_idx=v_temp==0;

    %include hets for now
    alt_idx=logical((v_temp==1)+(v_temp==-1));

    v_maf(i)=min([sum(ref_idx),sum(alt_idx)])/length(v_niche);

    %ref is major
    if sum(ref_idx)>=sum(alt_idx)

        v_niche_major=v_niche_code(ref_idx);
        v_niche_minor=v_niche_code(alt_idx);

        modal_niche_major=mode(v_niche_major);
        modal_niche_minor=mode(v_niche_minor);

        f_modal_niche_major(i)=sum(v_niche_major==modal_niche_major)/length(v_niche_major);
        f_modal_niche_minor(i)=sum(v_niche_minor==modal_niche_minor)/length(v_niche_minor);

        f_ferm_major(i)=sum(v_niche_major==3)/length(v_niche_major);
        f_ferm_minor(i)=sum(v_niche_minor==3)/length(v_niche_minor);

    elseif sum(ref_idx)<sum(alt_idx)

        v_niche_major=v_niche_code(alt_idx);
        v_niche_minor=v_niche_code(ref_idx);

        modal_niche_major=mode(v_niche_major);
        modal_niche_minor=mode(v_niche_minor);

        f_modal_niche_major(i)=sum(v_niche_major==modal_niche_major)/length(v_niche_major);
        f_modal_niche_minor(i)=sum(v_niche_minor==modal_niche_minor)/length(v_niche_minor);

        f_ferm_major(i)=sum(v_niche_major==3)/length(v_niche_major);
        f_ferm_minor(i)=sum(v_niche_minor==3)/length(v_niche_minor);

    end
    
    %also make matrix of niches to export
    for j=1:length(niche_names)
        
        niche_mat_ref(i,j)=sum(v_niche_code(ref_idx)==j);
        niche_mat_alt(i,j)=sum(v_niche_code(alt_idx)==j);
        
    end
    

end


chr_common=chr;
pos_common=pos;
alt_common=alt;
af_common=v_maf;

for i=1:length(alt_common)
    
    temp_str=strsplit(alt_common{i},',');
    alt_common{i}=temp_str{1};
    
end

%now need to match to structure stats for each mutation
load([dependency_directory '1K_data_annotated.mat'])

%filter on af
af(af>0.5)=af(af>0.5)-0.5;
af_idx=af>=0.05;

chr=chr(af_idx);
pos=pos(af_idx);
gene=gene(af_idx);
proteinEncoded=proteinEncoded(af_idx);

chr(cellfun(@isempty,chr))={'NA'};
for i=1:length(chr_common)
    
    if mod(i,10000)==0
        i
    end
    
    temp_chr_idx=ismember(chr,chr_common{i});
    
    temp_pos=pos(temp_chr_idx);
    temp_alt=alt(temp_chr_idx);
    temp_gene=gene(temp_chr_idx);
    temp_protein_encoded=proteinEncoded(temp_chr_idx);
    
    temp_pos_idx=temp_pos==pos_common(i);
    
    if sum(temp_pos_idx)>0
        
        temp_alt=temp_alt(temp_pos_idx);
        for j=1:length(temp_alt)
    
            temp_str=strsplit(temp_alt{j},',');
            temp_alt{j}=temp_str{1};

        end
        temp_gene=temp_gene(temp_pos_idx);
        temp_protein_encoded=temp_protein_encoded(temp_pos_idx);

        temp_alt_idx=ismember(temp_alt,alt_common{i});

        if sum(temp_alt_idx)>0
            
            gene_common{i}=temp_gene{temp_alt_idx};
            protein_encoded_common{i}=temp_protein_encoded{temp_alt_idx};
            
        end
        
    end
    
end


for i=1:length(protein_encoded_common)
    
    if ~isempty(protein_encoded_common{i})
        
        temp_str=strsplit(protein_encoded_common{i},'.');
        
        ref_res=temp_str{2}(1:3);
        alt_res=temp_str{2}((end-2):end);
        
        if ~strcmp(ref_res,alt_res) %missense
        
            residue_common(i)=str2num(temp_str{2}(4:(end-3)));
            
        end
        
    end
    
end

output_idx=residue_common>0;

chr_to_output=chr_common(output_idx);
pos_to_output=pos_common(output_idx);
gene_to_output=gene_common(output_idx);
residue_to_output=residue_common(output_idx);

af_to_output=af_common(output_idx);
niche_mat_ref_to_output=niche_mat_ref(output_idx,:);
niche_mat_alt_to_output=niche_mat_alt(output_idx,:);

%calculate p values here to save time


%calculate niche enrichment p values
niche_p_mat=nan(length(chr_to_output),length(niche_names));
niche_alt_enrichment_mat=niche_p_mat;
for i=1:length(chr_to_output)
    
    for j=1:length(niche_names)
    
        n1=niche_mat_ref_to_output(i,j);
        n2=sum(niche_mat_ref_to_output(i,1:(j-1)))+sum(niche_mat_ref_to_output(i,(j+1):end));
        n3=niche_mat_alt_to_output(i,j);
        n4=sum(niche_mat_alt_to_output(i,1:(j-1)))+sum(niche_mat_alt_to_output(i,(j+1):end));
        temp_table=table([n1;n3],[n2;n4],...
            'VariableNames',{'this niche','all other'},'RowNames',{'ref','alt'});
        [h,p,stats]=fishertest(temp_table);
        
        niche_p_mat(i,j)=p;
        niche_alt_enrichment_mat(i,j)=(n3/n1)/(n4/n2);
        
    end
    
    %output best p and enrichment for that SNP
    [min_p,min_idx]=min(niche_p_mat(i,:));
    
    niche_min_p(i)=min_p;
    niche_min_enrichment(i)=niche_alt_enrichment_mat(i,min_idx);
    niche_min_identitiy(i)=min_idx;
    
end


% histogram(-log10(length(chr_to_output)*niche_min_p))
% 
% gene_to_output(-log10(length(chr_to_output)*niche_min_p)>10)'
% 
% scatter(reshape(niche_alt_enrichment_mat,[],1),...
%     reshape(-log10(niche_p_mat),[],1))




%for table
%Xchr
%Xpos
%Xgene
%Xresidue #
%structure
%asa
%neighbors
%Xmaf
%niche p
%Xref count for each niche
%Xalt count for each niche

%look up structure, asa and neighbors
load([dependency_directory 'secondary_structure_data.mat'])
load([dependency_directory 'asa_data.mat'])
load([dependency_directory 'neighbor_data.mat'])


asa_to_output=nan(length(gene_to_output),1);
neighbors_to_output=nan(length(gene_to_output),1);
structure_to_output=nan(length(gene_to_output),1);
for i=1:length(gene_to_output)
    
    gene_idx=ismember(genes_to_use,gene_to_output{i});
    
    if sum(gene_idx)>0
        
        asa_to_output(i)=asa_mat(gene_idx,residue_to_output(i));
        neighbors_to_output(i)=neighbor_mat(gene_idx,residue_to_output(i));
        
        [~,temp_structure]=structure_types(secondary_mat{gene_idx,residue_to_output(i)});
        structure_to_output(i)=temp_structure;
        
    end
    
end


to_output=table(chr_to_output',pos_to_output',gene_to_output',...
    residue_to_output',structure_to_output,asa_to_output,neighbors_to_output,...
    af_to_output',niche_min_enrichment',-log10(niche_min_p'),niche_min_identitiy',niche_mat_ref_to_output,niche_mat_alt_to_output,...
    'VariableName',{'chr','pos','gene','residue','secondary','asa','neighbors',...
    'maf','niche_enrichment','niche_p_value','niche_id','ref_niches','alt_niches'});
writetable(to_output,[dependency_directory '1K_common_annotated_niche.csv'])



% 
% maf_thresh=0.2;
% p_thresh=10;
% 
% table_to_use=to_output(to_output.maf>maf_thresh,:);
% 
% sum(table_to_use.niche_p_value>p_thresh)
% 
% temp_idx=table_to_use.niche_p_value>p_thresh;
% 
% to_plot{1}=table_to_use.asa(temp_idx);
% to_plot{2}=table_to_use.asa(~temp_idx);
% 
% subplot(2,4,1)
% easy_box(to_plot)
% 
% 
% to_plot{1}=table_to_use.neighbors(temp_idx);
% to_plot{2}=table_to_use.neighbors(~temp_idx);
% 
% 
% subplot(2,4,2)
% easy_box(to_plot)





toc





