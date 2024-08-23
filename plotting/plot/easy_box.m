function []=easyBox(inputArray)

maxLength=max(cellfun(@length,inputArray));

matToPlot=nan(maxLength,length(inputArray));

for i=1:length(inputArray)
    
    matToPlot(1:length(inputArray{i}),i)=inputArray{i};
    
end

boxplot(matToPlot,'Symbol','.','Colors','k')

end