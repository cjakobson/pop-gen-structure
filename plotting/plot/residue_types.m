%function that returns structure type distribution
function [dist_out,array_out]=residue_types(inputArray)

aa_labels={'A','V','I','L','M','F','Y','W','S','T','N','Q','H','R','K','D','E','C','G','P'};


%allow 2D input
array_out=zeros(size(inputArray));

for i=1:length(aa_labels)

    array_out(ismember(inputArray,aa_labels{i}))=i;

end

for i=1:length(aa_labels)
     
    dist_out(i)=sum(sum(array_out==i));
    
end

dist_out=dist_out./sum(dist_out);

