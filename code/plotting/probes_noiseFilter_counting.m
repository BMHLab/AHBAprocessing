function figure3A
clear all; 
close all; 
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

signalThresholds = 0:0.01:1;
% count how many genes/probes are there and how many are excluded with
% signal threshold of 0.5
probes = DataTableProbe.ProbeName{1}; 
genes = DataTableProbe.EntrezID{1}; 

[uGenes, uGenesIND] = unique(genes);
fprintf('%d unique probes\n', length(probes))
fprintf('%d unique genes\n', length(uGenes))

%% find repeating entrezIDs and calculate variances for them, then take a probe of max variance.

% % ------------------------------------------------------------------------------
% First, find probes that have very noisy data and remove them from consideration
% Threshold for removing those proges is defined as the percentage of
% samples a probe has expression higher that background
% % ------------------------------------------------------------------------------
nrProbes = zeros(length(signalThresholds),1); 
nrGenes = zeros(length(signalThresholds),1); 

for j=1:length(signalThresholds)
    
signalThreshold = signalThresholds(j);     
    
signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=signalThreshold);
indRemoveProbes = find(signalLevel<signalThreshold);
keepProbes = signalLevel>=signalThreshold;

% remove selected probes from data and perform other calculations only on
% non-noisy probes
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);

nrProbes(j) = length(ProbeName); 
nrGenes(j) = length(unique(EntrezID)); 
end

fig=figure; 
set(fig,'defaultAxesColorOrder',[.92 .35 .24; 0 .5 .5]);

yyaxis left
plot(signalThresholds, nrGenes, '-o', ...
    'Color', [.6 .6 .6], 'LineWidth',1,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[.6 .6 .6],...
    'MarkerFaceColor',[.96 .63 .55]);
xlabel('Intensity-based filtering threshold'); 
ylabel('Number of genes')
yticklabels([0 5000 10000 15000 20000 25000])
yticks([0 5000 10000 15000 20000 25000])
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
set(gcf,'color','w');
hold on; 

yyaxis right
plot(signalThresholds, nrProbes, '-o', ...
    'Color', [.6 .6 .6], 'LineWidth',1,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[.6 .6 .6],...
    'MarkerFaceColor',[.59 .87 .82]);
ylabel('Number of probes','FontSize', 14)
yticklabels([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
yticks([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
%legend({'Number of genes', 'Number of probes'}); 
set(gca,'FontSize', 16)
box off
h = legend({'Number of genes', 'Number of probes'}); 
set(h,'fontsize',20)

