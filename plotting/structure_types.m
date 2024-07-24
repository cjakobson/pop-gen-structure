%function that returns variant type distribution
function [distOut,arrayOut]=structureTypes(inputArray)

%helices
variantArray{1}={'H'};
variantArray{2}={'G'};
variantArray{3}={'I'};
%sheets
variantArray{4}={'B'};
variantArray{5}={'E'};
%turns
variantArray{6}={'S'};
variantArray{7}={'T'};
%none
variantArray{8}={'U'};

arrayOut=zeros(1,length(inputArray));

for i=1:length(variantArray)

    arrayOut(ismember(inputArray,variantArray{i}))=i;

end

for i=1:length(variantArray)
     
    distOut(i)=sum(arrayOut==i);
    
end

distOut=distOut./sum(distOut);

