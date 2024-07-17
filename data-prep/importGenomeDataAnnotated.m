%scan text file and convert line-by-line
function []=importGenomeDataAnnotated(nRows)
tic
fID=fopen('/Volumes/cmjBackup/1011_matrix_roman_ann.vcf');
%fID=fopen('/Volumes/GeldocHDD/1011genomes/1011Matrix.ann.gvcf');


%N=ceil(10^logRows);
N=nRows;

m=1;

%this number is from snpEff
nMissenseSynonymous=585540+540294;

chr=cell(nMissenseSynonymous,1);
pos=nan(nMissenseSynonymous,1);
af=nan(nMissenseSynonymous,1);
gene=cell(nMissenseSynonymous,1);
type=cell(nMissenseSynonymous,1);
dnaEncoded=cell(nMissenseSynonymous,1);
proteinEncoded=cell(nMissenseSynonymous,1);


%need to implement check for end of file

for k=1:N
    
    txtLine=fgetl(fID);
    if ischar(txtLine)
        if ~strcmp(txtLine(1),'#')
            temp=strsplit(txtLine,'\t');
            if length(temp)>=10
                tempChr=temp{1};
                tempPos=str2double(temp{2});
                temp2=strsplit(temp{8},';');
                temp3=strsplit(temp2{2},{'=',','});
                tempAf=str2double(temp3{2});
                temp4=temp2{end};
                temp5=strsplit(temp4,'|');
                tempType=temp5{2};
                if strcmp(tempType,'missense_variant')||...
                        strcmp(tempType,'synonymous_variant')
                    tempGene=temp5{5};
                    tempDnaEncoded=temp5{10};
                    tempProteinEncoded=temp5{11};

                    chr{m}=tempChr;
                    pos(m)=tempPos;
                    af(m)=tempAf;
                    gene{m}=tempGene;
                    type{m}=tempType;
                    dnaEncoded{m}=tempDnaEncoded;
                    proteinEncoded{m}=tempProteinEncoded;
                    m=m+1;

                    if mod(m,1000)==0
                        m
                        k
                    end


                end
                %m=m+1;
            end
        end
    end
end
fclose(fID);

% histogram(log10(af(ismember(type,'missense_variant'))),-4:0.4:0)
% hold on
% histogram(log10(af(ismember(type,'synonymous_variant'))),-4:0.4:0)

save('1002dataAnn.mat')

toc
end






