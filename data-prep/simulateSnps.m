function outputTable=simulateSnps(systematicName,dnaSequence)

%some sequences seem to not be mod(length,3)==0
nBasesToUse=3*floor(length(dnaSequence)/3);

baseArray={'A','C','G','T'};

initialBase=cell(3*nBasesToUse,1);
newBase=initialBase;
initialResidue=initialBase;
newResidue=initialBase;
initialCodon=initialBase;
newCodon=initialBase;

isMis=zeros(3*nBasesToUse,1);
residueNumber=isMis;

isTs=isMis;

m=1;
for i=1:nBasesToUse

    tempBase=dnaSequence(i);
    residueIdx=(1:3)+(ceil(i/3)-1)*3;
    tempResidue=nt2aa(dnaSequence(residueIdx),'AlternativeStartCodons',false);
    
    mutationArray=baseArray(~ismember(baseArray,tempBase));
    for j=1:length(mutationArray)
        
        initialBase{m}=tempBase;
        newBase{m}=mutationArray{j};
        
        newSequence=dnaSequence;
        newSequence(i)=mutationArray{j};
        
        initialResidue{m}=tempResidue;
        newResidue{m}=nt2aa(newSequence(residueIdx),'AlternativeStartCodons',false);
        
        initialCodon{m}=dnaSequence(residueIdx);
        newCodon{m}=newSequence(residueIdx);
        
        isMis(m)=~strcmp(initialResidue{m},newResidue{m});
        
        %Ts/Tv
        if (strcmp(initialBase{m},'A')&&strcmp(newBase{m},'G'))||...
                (strcmp(initialBase{m},'G')&&strcmp(newBase{m},'A'))||...
                (strcmp(initialBase{m},'C')&&strcmp(newBase{m},'T'))||...
                (strcmp(initialBase{m},'T')&&strcmp(newBase{m},'C'))
            isTs(m)=1;
        end
        
        residueNumber(m)=ceil(i/3);
        
        m=m+1;
        
    end
    
end

outputTable=table(residueNumber,initialBase,newBase,initialResidue,newResidue,...
    initialCodon,newCodon,isMis,isTs);

save(['mutationOutput/' systematicName 'mutationTable.mat'],'outputTable')

end