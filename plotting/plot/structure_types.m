%function that returns structure type distribution
function [dist_out,array_out]=structure_types(inputArray)

%helices
variant_array{1}={'H'};
variant_array{2}={'G'};
variant_array{3}={'I'};
%sheets
variant_array{4}={'B'};
variant_array{5}={'E'};
%turns
variant_array{6}={'S'};
variant_array{7}={'T'};
%none
variant_array{8}={'U'};

%allow 2D input
array_out=zeros(size(inputArray));

for i=1:length(variant_array)

    array_out(ismember(inputArray,variant_array{i}))=i;

end

for i=1:length(variant_array)
     
    dist_out(i)=sum(sum(array_out==i));
    
end

dist_out=dist_out./sum(dist_out);

