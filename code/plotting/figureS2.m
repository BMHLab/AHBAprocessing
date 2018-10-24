function [rDS, rRNA] = figureS3()
% plot average corr to RNA-seq vs signal proportion
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat'); 
% get noise values and sort them to match the order of avgCorr
numBins = 1000; 
noiseVal = probeInformationALL.signalProp; 
corrVal = cell2mat(avgCorr);
bins = linspace(0,1,numBins+1); 

for i=1:length(bins)-1 
    ind = find(noiseVal>bins(i) & noiseVal<=bins(i+1)); 
    meanCorr(i) = nanmean(corrVal(ind)); 
end
[rRNA,pRNA] = corr(noiseVal, corrVal, 'rows', 'complete','type', 'Spearman'); 
 
figure; set(gcf,'Position',[300 300 1300 500])
f1 = subplot(1,2,1); 
plot(bins(1:numBins), meanCorr, 'o', 'MarkerSize', 10 ,'MarkerFaceColor', [.7 .7 .7], 'MarkerEdgeColor', [.5 .5 .5]); 
xlabel({'The proportion of samples with', 'expression exceeding the background'})
ylabel({'Average correlation between', 'microarray and RNAseq'})
set(gca,'FontSize', 16)
set(gcf,'color','w'); 
box(f1,'off')
 
% plot DS vs signal proportion
load('100DS82scaledRobustSigmoidNSGRNAseqnoQC1Lcortex_ROI_NOdistCorrEuclidean.mat'); 
DS = probeInformation.DS; 
probeID = probeInformation.ProbeID;  

% for each DS threshold get avCorr
for g=1:length(probeID)
    indNOISE = find(probeInformationALL.ProbeID==probeID(g)); 
    indDS = find(probeID==probeID(g)); 
    mxNoise(g)  = probeInformationALL.signalProp(indNOISE); 
    DSval(g) = DS(indDS); 
end

for i=1:length(bins)-1 
    ind = find(mxNoise>bins(i) & mxNoise<=bins(i+1)); 
    DSmean(i) = nanmean(DSval(ind)); 
end

[rDS,pDS] = corr(mxNoise', DSval', 'rows', 'complete','type', 'Spearman'); 
f2 = subplot(1,2,2); 
plot(bins(1:numBins),  DSmean, 'o', 'MarkerSize', 10 ,'MarkerFaceColor', [.7 .7 .7], 'MarkerEdgeColor', [.5 .5 .5]); 
xlabel({'The proportion of samples with', 'expression exceeding the background'})
ylabel({'Average DS'})
set(gca,'FontSize', 16)
set(gcf,'color','w'); 
box(f2,'off')
end
