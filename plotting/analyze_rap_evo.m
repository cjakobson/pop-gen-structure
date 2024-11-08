
clear

tic

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

figure_counter=1;

code_directory=[filebase 'Documents/GitHub/pop-gen-structure/'];
dependency_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/pop-gen-structure-output/'];
addpath([code_directory 'data-prep'])
addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/plot'])

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



%dotplot of final adapted growth rates from pilot
strains={'BY4741','BY4741Dmsh6','BY4743'};%,'BY4743DDmsh6'};

clear to_plot
for i=1:length(strains)

    spot_size_table{i}=readtable([dependency_directory 'gitter-data/' strains{i} '_pilot.csv']);
    spot_size_mat{i}=table2array(spot_size_table{i}(:,2:end));

    to_plot{i}=spot_size_mat{i}(end,:);

end




figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,4,1)
hold on
easy_dotplot(to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(strains)
ylabel('final spot size')
axis square
title('rapamycin pilot')

evo_thresh_spot_size=1500;
plot(xlim,[evo_thresh_spot_size evo_thresh_spot_size],':r')


clear n_to_plot
for i=1:length(to_plot)

    n_to_plot(i)=sum(to_plot{i}>evo_thresh_spot_size);
    
end
subplot(2,8,3)
hold on
bar(n_to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(strains)
ylabel('number of adapted lineages')
%axis square


%median normalize to non-adapted
for i=1:length(to_plot)

    to_plot{i}=to_plot{i}/median(to_plot{i},'omitnan');

end


subplot(2,4,3)
hold on
easy_dotplot(to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(strains)
ylabel('final spot size (norm.)')
axis square
title('rapamycin pilot')

evo_thresh_spot_size_norm=2;
plot(xlim,[evo_thresh_spot_size_norm evo_thresh_spot_size_norm],':r')


clear n_to_plot
for i=1:length(to_plot)

    n_to_plot(i)=sum(to_plot{i}>evo_thresh_spot_size_norm);
    
end
subplot(2,8,7)
hold on
bar(n_to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(strains)
ylabel('number of adapted lineages')
%axis square




set(gcf,'PaperPositionMode','auto')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



%plot example curves over time
condition_names={'YPD+rap1','YPD+rap2','YPDonly'};
start_passage=2;

for i=1:length(condition_names)
    
    condition_spot_size_table{i}=readtable([dependency_directory 'gitter-data/' condition_names{i} '_evo.csv']);
    condition_spot_size_mat{i}=table2array(condition_spot_size_table{i}(:,2:end));



    condition_slope_table{i}=readtable([dependency_directory 'gitter-data/' condition_names{i} '_evo_slope.csv']);
    condition_slope_mat{i}=table2array(condition_slope_table{i});

    condition_slope_mat{i}=condition_slope_mat{i}(start_passage:end,:);


end



figure('units','normalized','outerposition',[0 0 1 1])

evo_thresh_slope1=1;
evo_thresh_slope2=0.5;


for i=1:2%length(condition_names)
    
    

    [n_time_points,n_strains]=size(condition_slope_mat{i});

    subplot(2,2,i)
    hold on

    v_final=condition_slope_mat{i}(end,:);
    adapted_idx{i}=find(v_final>evo_thresh_slope1);
    non_adapted_idx{i}=find(v_final<=evo_thresh_slope2);

    %example adapted and non-adapted trajectories

    for j=1:length(adapted_idx{i})

        v1=1:n_time_points;
        v2=condition_slope_mat{i}(:,adapted_idx{i}(j));

        plot(v1,v2,'-r')
        scatter(v1,v2,20,'r','filled')

    end

    for j=101:200%length(adapted_idx)

        v1=1:n_time_points;
        v2=condition_slope_mat{i}(:,non_adapted_idx{i}(j));

        plot(v1,v2,'-k')
        scatter(v1,v2,20,'k','filled')

    end

    ylim([0 2.5])
    xlim([0.5 n_time_points+0.5])
    %axis square
    title(condition_names{i})

end
ylabel('growth rate')
xlabel('passage #')




%skip first passage -- unusually high growth due to no drug before
start_time=24;
for i=1:2%length(condition_names)

    [n_time_points,n_strains]=size(condition_spot_size_mat{i});

    subplot(2,2,i+2)
    hold on

    v1=1:n_time_points;
    v2=condition_spot_size_mat{i}(:,non_adapted_idx{i}(1));

    plot(v1(start_time:end),v2(start_time:end),'-k')
    scatter(v1(start_time:end),v2(start_time:end),20,'k','filled')



    v1=1:n_time_points;
    v2=condition_spot_size_mat{i}(:,adapted_idx{i}(11));

    plot(v1(start_time:end),v2(start_time:end),'-r')
    scatter(v1(start_time:end),v2(start_time:end),20,'r','filled')
    
    ylim([0 3000])
    title(condition_names{i})
    ylabel('spot size')
    xlabel('time (hrs)')


end





set(gcf,'PaperPositionMode','auto')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;



figure('units','normalized','outerposition',[0 0 1 1])

clear to_plot
%dot plots of final slopes
for i=1:2%length(condition_names)
    
    to_plot{i}=condition_slope_mat{i}(end,:);

    n_adapted(i)=sum(to_plot{i}>evo_thresh_slope1);

end

subplot(2,8,1)
hold on
easy_dotplot(to_plot)
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(condition_names(1:length(to_plot)))
ylabel('final passage growth rate')
%axis square
title('rapamycin evolution')
ylim([0 2.5])



subplot(2,8,2)
hold on
bar(n_adapted)
ylim([0 100])
ylabel('# adapted clones')
xticks(1:length(to_plot))
xtickangle(45)
xticklabels(condition_names(1:length(to_plot)))



set(gcf,'PaperPositionMode','auto')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;







%compare final slope with re-phenotyping of slope in rap and YPD only
exp_name={'round1','round2','round3','round4'};

plate_labels={'rap1','rap2','YPDonly'};

condition_names={'rap','YPD'};


for i=1:length(exp_name)
    
    for j=1:length(plate_labels)
        
        if i<=2
            k=1;
        else
            k=2;
        end
    

        phenotyping_slope_mat{i,j}=table2array(readtable([dependency_directory 'gitter-data/'...
            exp_name{i} '_' plate_labels{j} '=>' condition_names{k} '_phen_slope.csv']));

        phenotyping_slope_names{i,j}=[exp_name{i} '_' plate_labels{j} '=>' condition_names{k}];

    end

end



figure('units','normalized','outerposition',[0 0 1 1])


%scatter replicates
for i=1:length(plate_labels)

    v1=phenotyping_slope_mat{1,i};
    v2=phenotyping_slope_mat{2,i};

    phenotyping_mean_rap(i,:)=mean([v1; v2]);

    subplot(2,4,i)
    hold on
    scatter(v1,v2,10,'k','filled')
    title([plate_labels{i} ' in ' condition_names{1}])
    xlim([0 2])
    ylim(xlim)
    axis square
    xlabel('replicate 1')
    ylabel('replicate 2')

end

for i=1:length(plate_labels)

    v1=phenotyping_slope_mat{3,i};
    v2=phenotyping_slope_mat{4,i};

    phenotyping_mean_ypd(i,:)=mean([v1; v2]);

    subplot(2,4,i+4)
    hold on
    scatter(v1,v2,10,'k','filled')
    title([plate_labels{i} ' in ' condition_names{2}])
    xlim([0 2])
    ylim(xlim)
    axis square
    xlabel('replicate 1')
    ylabel('replicate 2')

end




set(gcf,'PaperPositionMode','auto')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;






%no cross-adaptation
figure('units','normalized','outerposition',[0 0 1 1])

for i=1:length(plate_labels)

    v1=phenotyping_mean_rap(i,:);
    v2=phenotyping_mean_ypd(i,:);

    subplot(2,4,i)
    hold on
    scatter(v1,v2,10,'k','filled')
    title(plate_labels{i})
    xlim([0 2])
    ylim(xlim)
    axis square
    xlabel('growth rate in rap')
    ylabel('growth rate in ypd')

end


adapted_slope_thresh=0.75;


temp_labels={'non-adapted in rap',...
    'adapted in rap',...
    'non-adapted in YPD',...
    'adapted in YPD'};
for i=1:length(plate_labels)

    v1=phenotyping_mean_rap(i,:);
    v2=phenotyping_mean_ypd(i,:);

    adapted_idx=v1>adapted_slope_thresh;

    clear to_plot

    to_plot{1}=v1(~adapted_idx);
    to_plot{2}=v1(adapted_idx);

    to_plot{3}=v2(~adapted_idx);
    to_plot{4}=v2(adapted_idx);

    subplot(2,4,i+4)
    easy_box(to_plot)
    title(plate_labels{i})

    xticks(1:length(to_plot))
    xtickangle(45)
    xticklabels(temp_labels)
    ylabel('growth rate')
    ylim([0 2])
    axis square
    for j=1:length(to_plot)
        text(j,0.1,num2str(length(to_plot{j})))
    end

end




set(gcf,'PaperPositionMode','auto')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-dsvg','-r0')
print([output_directory 'rap_evo_figure_' num2str(figure_counter)],'-djpeg','-r300')
figure_counter=figure_counter+1;





toc




