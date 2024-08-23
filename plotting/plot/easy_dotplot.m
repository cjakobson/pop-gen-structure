function []=easy_dotplot(input_array)

hold on
for i=1:length(input_array)

    v1=i+0.1*randn(length(input_array{i}),1);
    v2=input_array{i};

    scatter(v1,v2,10,'k','filled','MarkerFaceAlpha',0.5)

end

end