function outputTable=getNeighbors(systematicName,pdbId)

fileToGet=['/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/UP000002311_559292_YEAST/' pdbId '_10A.txt'];
    
if exist(fileToGet)

    %get neighbors within 10Ang
    tempNeighbors=readtable(fileToGet);
    
    outputTable=tempNeighbors;

    save(['mutationOutput/' systematicName 'neighborTable.mat'],'outputTable')
    
end

    





end