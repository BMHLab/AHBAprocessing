% correlate variance and intensity for each probe - we expect negative
% correlation
function figureS1()

cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
figure; set(gcf,'Position',[300 300 1000 500])
set(gcf,'color','w');
% calculate;
for j=1:2
    if j==1
        whatMeasure='variance';
        nameLabel = 'variance';
        numBins = 500;
    elseif j==2
        whatMeasure='CV';
        nameLabel = 'coeficient of variation';
        numBins = 250;
    end
    
    for i=1:size(Expressionall,1)
        if strcmp(whatMeasure, 'CV')
            V(i) = (std(Expressionall(i,:)))/mean(Expressionall(i,:));
        elseif strcmp(whatMeasure, 'variance')
            V(i) = var(Expressionall(i,:));
        end
        I(i) = mean(Expressionall(i,:));
        
    end
    subplot(1,2,j)
    [xThresholds,yMeans] = BF_PlotQuantiles(I,V,numBins,0,0);
    xlabel('Mean probe intensity');ylabel(sprintf('Probe %s', nameLabel));
    set(gca,'FontSize', 18)
   
end
    cd ../../..
end
