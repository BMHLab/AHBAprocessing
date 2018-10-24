function figure4B(QClabel)
% script to evaluate correlations between RNAseq and microarray data and
% make list for enrichment
cd ('data/genes/processedData')
load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq%s.mat', QClabel))
    
ind = cell2mat(cellfun(@(x)any(isnan(x)),avgCorr,'UniformOutput',false)); 
fprintf('%d genes are overlapping between RNA-seq and microarray datasets\n', length(find(ind==0)))

avgCorr = avgCorr(ind==0); 

for i=1:length(find(ind==0))
    maxCor(i) = max(avgCorr{i}); 
end

maxCor = maxCor'; 
% get min and max values
minimal = min(maxCor); 
maximal = max(maxCor); 
med = median(maxCor); 
% how many <0.3
low = length(find(maxCor<0.3))/(length(maxCor));
% how many >0.5
high = length(find(maxCor>0.5))/(length(maxCor));


figure; colors = [.96 .63 .55; 1 .46 .22]; 
histogram(maxCor, 100,'EdgeColor',[.6 .6 .6],...
    'FaceColor',colors(1,:)); 
xtitleText = {'Average correlation between microarray and RNA-seq', 'gene expression'};
xlabel(xtitleText)
ylabel('Number of genes')
set(gcf,'color','w'); hold on; 

set(gca,'FontSize', 16)
box off

legendText = sprintf('%d genes', length(find(ind==0))); 
h = legend(legendText); 
set(h,'fontsize',20)
cd ../../..
end


