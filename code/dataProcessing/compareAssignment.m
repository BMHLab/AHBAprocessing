% import file as cell called dataCell
% create empty variables to store info

dataCell = importProbesALL('probes2annotateALL_merged_readAnnotation.txt');
probeOK = zeros(size(dataCell,1)-1,1); 
geneNames = cell(size(dataCell,1)-1,1); 
probeNames = cell(size(dataCell,1)-1,1);
% for each row select non empty elements
for i=1:length(probeOK)
   % take part where assignment was done tratring from 6th column
   probe = dataCell(i,6:end);
   % count how many matches there are
   indCell = dataCell(i,5); 
   numMatches = length(strfind(probe{1}, ';')); 
   %numMatches = indCell{1}; 
   % if there is one match, say that gene probe is OK and save gene name by
   % splitting string to select relevant part
   if numMatches==1
       probeOK(i) = 1;
       C = strsplit(probe{1},','); 
       B = strsplit(C{2},'|'); 
       geneNames{i} = B{1};
       probeNames{i} = dataCell{i}; 
   else
   % for each nonempty probe split it into two parts based on comma and
   % check if gene name is the same for all matches
   geneName = cell(numMatches,1);
   C = strsplit(probe{1},',');
   for p=1:numMatches
          B = strsplit(C{p+1},'|'); 
          geneName{p} = B{1};
   end
     % check if all match
   doMatch = zeros(numMatches,numMatches); 
   for j=1:numMatches-1
       for k=j+1:numMatches
           doMatch(j,k) = strcmp(geneName{j}, geneName{k});
       end
   end
   % check how many match
   numMatchesName = sum(doMatch(:)); 
   numMax = (numMatches.*numMatches - numMatches)/2; 
   % if all gene names match - good to go
   probeOK(i) = numMatchesName==numMax; 
   if probeOK(i)==1
       geneNames{i} = B{1}; 
       probeNames{i} = dataCell{i}; 
   end
   end
   
end

% remove rows where probes were not mapped uniquely
probeNames(probeOK==0) = []; 
geneNames(probeOK==0) = [];
% combine gene names and probe names to one table
hg38match = table(probeNames,geneNames); 
% take data downloaded from
% ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/ about gene
% information (gene names mapped to IDs). 
% import first 4 columns into a table and for each gene in the hg38match
% table find a corresponding gene ID. 
geneData = importGeneFile('Homo_sapiens_gene_info.txt');
numGenes = size(hg38match,1); 

for g=1:numGenes
    gene = geneNames{g}; % gene names from converted probes

        x = strcmp(geneData(:,2), gene); 
        k = find(x); 
        if length(k)==1
            hg38match.ID(g) = geneData{k,1};
        else
            hg38match.ID(g) = NaN;
        end

end

allenProbes = importAllenProbes('Probes.xlsx'); 
% import original data from allen
%fprintf(1,'Removing CUST probes\n')
%cust = strfind(allenProbes.probeNames, 'CUST');
%remInd = find(~cellfun(@isempty,cust));
%fprintf(1,'%d CUST probes removed\n', length(remInd))
%allenProbes(remInd,:) = [];

% for all probes that were re-annotated check if entrezID matches the
% reannotated entrezID.
hg38match = comparehg38VSAllen(hg38match, allenProbes); 
 
fprintf('%d probes match\n', length(find(hg38match.compare==1))); 

fprintf('%d probes not assigned an ID in NCBI - removed\n', length(find(hg38match.compare==3))); 
%indREM = find(hg38match.compare==3); 
%hg38match(indREM,:) = []; 

fprintf('%d probes changed\n', length(find(hg38match.compare==2))); 
fprintf('%d probes do not match\n', length(find(hg38match.compare==0)));


fprintf('%d probes are mapped to genes\n', length(find(hg38match.compare==0))+ ...
    length(find(hg38match.compare==1))+length(find(hg38match.compare==2)))

save('reannotatedProbes.mat', 'hg38match')
%allenProbes = compareAllenVShg38(hg38match, allenProbes); 




