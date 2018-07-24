%% script to exclude probes that were assigned to multiple genes
% Aurina 2018-01-09


%% import allen probes and remove CUST probes as they can't be annotated
allenProbes = importAllenProbes('Probes.xlsx'); 
% import original data from allen
fprintf(1,'Removing CUST probes\n')
cust = strfind(allenProbes.probeNames, 'CUST');
remInd = find(~cellfun(@isempty,cust));
fprintf(1,'%d CUST probes removed\n', length(remInd))
allenProbes(remInd,:) = [];
%% make a txt file for all allen probes (to be imported to website)
% import a list to the website: http://asia.ensembl.org/biomart/martview/39c608cf8abba1c5248dbe73bcd9c639
% and save output as mart_export_updatedProbes.txt

%% load re-annotated probes
updatedProbes = importProbes('mart_export_updatedProbes.txt'); 

% find probes that are mentioned more than once
[~, i, j] = unique(updatedProbes.probeNames,'first');
indexToDupes = find(not(ismember(1:numel(updatedProbes.probeNames),i))); 

% exclude those probes 
updatedProbes(indexToDupes,:) = []; 
% now each probe is mapped to one gene
% It's time to compare allen probes with re-annotated probes. 
% keeping in mind that:
%1. probes that were not reannotated can't be verified
%2. probes that had no entrez IDs to begin with could be re-annotated and updated
% it makes moe sense to take all re-annotated probes and compare them to
% what allen supplied (than the other way around). 

updatedMatching = comparehg38VSAllen(updatedProbes, allenProbes); 

%% sumarise the information in numbers
nm = length(find(updatedMatching.compare==1));
fprintf(1,'%d probes are matching\n', nm)
nmm = length(find(updatedMatching.compare==0));
fprintf(1,'%d probes are mismatching\n', nmm)
nu = length(find(updatedMatching.compare==2));
fprintf(1,'%d probes are introduced with IDs\n', nu)

