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



%now need to match to structure stats for each mutation
%based on chr/pos/ref/alt



end


