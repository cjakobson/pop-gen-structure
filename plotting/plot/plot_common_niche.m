function [] = plot_common_niche(dependency_directory)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

ecotype_data=readtable([dependency_directory '1002_genomes_ecotypes.txt']);


load([dependency_directory '1002_data_common.mat'])


strain_names=ecotype_data.StandardID;

%match strain names
for i=1:length(strainString)

    temp_idx=find(ismember(strain_names,strainString{i}));

    if ~isempty(temp_idx)
            
        v_niche{i}=ecotype_data.Ecotype{temp_idx};

    else

        v_niche{i}='NA';

    end

end
niche_names=unique(v_niche);

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




figure('units','normalized','outerposition',[0 0 1 1])

for i=1:length(niche_names)
    v_to_plot(i)=sum(v_niche_code==i);
end
subplot(2,4,1)
bar(v_to_plot)
xticks(1:length(niche_names))
xtickangle(45)
xticklabels(niche_names)
axis square
ylabel('frequency')
title('broad niches')




v_bins=0.05:0.05:0.5;
v1=v_maf;
v2=f_ferm_major;
clear to_plot
for i=1:(length(v_bins)-1)
    temp_idx=logical((v1>=v_bins(i)).*(v1<v_bins(i+1)));
    to_plot{i}=v2(temp_idx);
end

subplot(2,4,2)
%scatter(v_maf,f_modal_niche_major,1,'k','filled')
easy_box(to_plot)
axis square
xticklabels(v_bins)
title('strains with major allele')
ylabel('f_{ferm}')
xlabel('MAF')
ylim([0 1])



v1=v_maf;
v2=f_ferm_minor;
clear to_plot
for i=1:(length(v_bins)-1)
    temp_idx=logical((v1>=v_bins(i)).*(v1<v_bins(i+1)));
    to_plot{i}=v2(temp_idx);
end

subplot(2,4,3)
%scatter(v_maf,f_modal_niche_major,1,'k','filled')
easy_box(to_plot)
axis square
xticklabels(v_bins)
title('strains with minor allele')
ylabel('f_{ferm}')
xlabel('MAF')
ylim([0 1])




chr_common=chr;
pos_common=pos;
alt_common=alt;

for i=1:length(alt_common)
    
    temp_str=strsplit(alt_common{i},',');
    alt_common{i}=temp_str{1};
    
end

%now need to match to structure stats for each mutation
load([dependency_directory '1K_data_annotated.mat'])

%filter on af
af(af>0.5)=af(af>0.5)-0.5;
af_idx=af>0.05;

chr=chr(af_idx);
pos=pos(af_idx);
gene=gene(af_idx);
proteinEncoded=proteinEncoded(af_idx);

chr(cellfun(@isempty,chr))={'NA'};
for i=1:length(chr_common)
    
    if mod(i,1000)==0
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


%for table
%chr
%pos
%gene
%residue #
%ref residue
%alt residue
%maf
%ref count for each niche
%alt count for each niche




end


