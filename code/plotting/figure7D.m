
function figure8D()

cd ('data/genes/processedData');
% load pre-computed distances on the surface and within GM volume
% cortical samples mapped to aparcaseg parcellation
load('DistancesONsurfaceXXX.mat'); 
surf = maskuHalf(distSamples);
surf(isnan(surf)) = []; 
load('distancesGM_MNIXXX.mat'); 
GM = maskuHalf(distSamples); 
GM(isnan(GM)) = []; 

% calculate Euclidean distances
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqQC82DistThresh2.mat')
ROIs = 1:34;
coordinates = vertcat(DataCoordinatesMNI{1}, DataCoordinatesMNI{2},DataCoordinatesMNI{3}, ...
    DataCoordinatesMNI{4},DataCoordinatesMNI{5},DataCoordinatesMNI{6}); 

[sampleIND,x] = find(coordinates(:,2)==ROIs); 
distSamples = pdist2(coordinates(sampleIND,3:end),coordinates(sampleIND,3:end));  
eucl = maskuHalf(distSamples); 
eucl(isnan(eucl)) = []; 

figure; 
set(gcf,'color','w'); 
histogram(surf, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor', [.75 .5 .51]); hold on;
histogram(GM, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor', [.85 .54 .4]); hold on; 
histogram(eucl, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor',[.43 .6 .47] ); hold on; 
legend({'On the surface', 'Within grey matter volume', 'Euclidean'})
box('off')
xlabel('Distance between samples (mm)'); 
ylabel('Number of sample pairs')
set(gca,'fontsize',15)
cd ../../..

end

