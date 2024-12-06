function []=easyBox(inputArray)

maxLength=max(cellfun(@length,inputArray));

matToPlot=nan(maxLength,length(inputArray));

for i=1:length(inputArray)
    
    matToPlot(1:length(inputArray{i}),i)=inputArray{i};
    
end

hold on
boxplot(matToPlot,'Symbol','','Colors','k')
for i=1:length(inputArray)
    
    v1=matToPlot(:,i);
    v2=0.05*randn(1,length(v1))+i;
    scatter(v2,v1,10,'k','filled')
    
end

end