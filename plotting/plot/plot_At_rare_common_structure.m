function [] = plot_At_rare_common_structure(dependency_directory)

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

[~,at_structure]=structure_types(at_data.secondary);

af_thresh=0.001;

temp_idx1=at_af<=af_thresh;

structure_labels={'alphahelix','310helix','pihelix','betasheet','betaladder',...
     'bend','turn','unstr.'};
 

for i=1:length(structure_labels)

    temp_idx2=at_structure==i;
    
    v1=at_asa(logical(temp_idx1.*temp_idx2));
    v2=at_asa(logical(~temp_idx1.*temp_idx2));
    
    to_plot1(i)=mean(v1,'omitnan')/mean(v2,'omitnan');

    [h p_val1(i)]=ttest2(v1,v2);

    
end



for i=1:length(structure_labels)

    temp_idx2=at_structure==i;
    
    v1=at_neighbor(logical(temp_idx1.*temp_idx2));
    v2=at_neighbor(logical(~temp_idx1.*temp_idx2));
    
    to_plot2(i)=mean(v1,'omitnan')/mean(v2,'omitnan');

    [h p_val2(i)]=ttest2(v1,v2);

    
end



hold on
% bar(to_plot,'BaseValue',1)
% for i=1:length(p_val)
%     text(i,0.98,num2str(p_val(i)),'Rotation',-45)
% end
bar([to_plot1; to_plot2]','BaseValue',1)
for i=1:length(p_val1)
    text(i,to_plot1(i),num2str(p_val1(i)),'Rotation',45)
    text(i,to_plot2(i),num2str(p_val2(i)),'Rotation',-45)
end
ylim([0.8 1.1])
%title('C_\alpha within 10 Ang.')
legend({'ASA (Ang.^2)','C_\alpha within 10 Ang.'})
xticks(1:length(structure_labels))
xtickangle(45)
xticklabels(structure_labels)
ylabel('rare/common')
axis square



end


