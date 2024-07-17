function outputTable=getDssp(systematicName,pdbId,dnaSequence)

fileToGet=['/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/UP000002311_559292_YEAST/AF-' pdbId '-F1-model_v1.pdbdssp.txt'];
    
if exist(fileToGet)

    %get DSSP
    tempDssp=readtable(fileToGet,'ReadVariableNames',0,'delimiter',' A ');

    secondary=cell(height(tempDssp),1);
    residue=cell(height(tempDssp),1);
    noStruct=zeros(height(tempDssp),1);
    sasa=zeros(height(tempDssp),1);
    phi=zeros(height(tempDssp),1);
    psi=zeros(height(tempDssp),1);
    posInOrf=zeros(height(tempDssp),1);
    for j=1:height(tempDssp)
        tempStr=tempDssp.Var2{j};
        if length(tempStr)>=4
            residue{j}=tempStr(1);
            secondary{j}=tempStr(4);
            if strcmp(tempStr(4),' ')
                noStruct(j)=1;
            end
            if length(tempStr)>=25
                sasa(j)=str2num(tempStr(23:25));
                phi(j)=str2num(tempStr(92:96));
                psi(j)=str2num(tempStr(98:102));
            end
        end
        posInOrf(j)=j/height(tempDssp);
    end

    %calculate lengths of runs in secondary structure (to
    %find positions in 'metahelix', etc
    tempIsRun=zeros(length(secondary),1);
    for j=2:length(secondary)
        if strcmp(secondary{j},secondary{j-1})
            tempIsRun(j-1)=1;
        end
    end

    tempPosInRun=zeros(length(secondary),1);
    tempCounter=0;
    for j=1:length(tempPosInRun)
        if tempIsRun(j)>0
            tempCounter=tempCounter+tempIsRun(j);
            tempPosInRun(j)=tempCounter;
        else
            tempCounter=0;
        end
    end

    tempLengthOfRun=tempPosInRun;
    for j=1:(length(tempLengthOfRun)-1)
        tempIdx=length(tempLengthOfRun)-(j-1);
        if tempPosInRun(tempIdx)==0
            tempLengthOfRun(tempIdx)=0;
        elseif tempLengthOfRun(tempIdx)>tempLengthOfRun(tempIdx-1)
            tempLengthOfRun(tempIdx-1)=tempLengthOfRun(tempIdx);
        end
    end
    
    for j=1:length(secondary)
        if strcmp(secondary{j},' ')
            secondary{j}='U';    %use U for unstructured (as opposed to blank which is missing data)
        end
    end

    posInRun=tempPosInRun;
    lengthOfRun=tempLengthOfRun;
    
    outputTable=table(residue,secondary,sasa,phi,psi,posInOrf,posInRun,lengthOfRun);

    save(['mutationOutput/' systematicName 'dsspTable.mat'],'outputTable','dnaSequence')
    
end

    





end