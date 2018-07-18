function figure6E(numROIs)
cd ('data/genes/processedData'); 
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqQC82DistThresh2.mat')

selectRegions = 1:numROIs;
D = cell(6,2);
subjNr = cell(6,1);

for s=1:6
data = DataExpression{s}; 
    select = ismember(data(:,2), selectRegions);
    cortexData = data(select==1,3:end);

        expData1 = BF_NormalizeMatrix(cortexData,'scaledRobustSigmoid');
        expData2 = BF_NormalizeMatrix(cortexData,'zscore');
         

    D{s,1} = expData1; 
    D{s,2} = expData2;
    subjNr{s} = data(select==1,1);
    R{s} = data(select==1,2);
 
end

expressionSRS = vertcat(D{1,1}, D{2,1}, D{3,1}, D{4,1}, D{5,1}, D{6,1}); 
expressionZscore = vertcat(D{1,2}, D{2,2}, D{3,2}, D{4,2}, D{5,2}, D{6,2}); 
subjects = vertcat(subjNr{1}, subjNr{2}, subjNr{3}, subjNr{4}, subjNr{5}, subjNr{6}); 

C = zeros(length(subjects),3); 
for s=1:length(subjects)
    if subjects(s)==1
        C(s,:) = [0 0 .5]; % dark blue
    elseif subjects(s)==2
        C(s,:) = [.53 .83 .97]; % light blue
    elseif subjects(s)==3
        C(s,:) = [.7 .7 .7]; % grey
    elseif subjects(s)==4
        C(s,:) = [1 .60 .40]; % orange
    elseif subjects(s)==5
        C(s,:) = [.81 .07 .15]; % red
    elseif subjects(s)==6
        C(s,:) = [.99 .87 .09]; % yellow
    end
end

sz = 60; 
geneIND = 5127; 


figure;set(gcf,'color','w');
scatter(expressionZscore(:,geneIND),expressionSRS(:,geneIND),sz,C,'filled','MarkerEdgeColor',[.55 .55 .55],'LineWidth',1);
xlabel('z-score', 'FontSize', 20); ylabel('SRS', 'FontSize', 20); 
set(gca,'FontSize', 18)
title(probeInformation.GeneSymbol{geneIND}, 'FontSize', 20); 
cd ../../..
end