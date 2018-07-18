% Calculate distances within GM
% 2018/03/29
% load the data you are interested in
parcellation = 82; 
distThreshold = 2; 
% this will load samples that are assigned to parcellation at a set
% threshold -each parcellation is slightly different, so the number of
% samples assigned is slightly different too. - need to do this for each
% parcellation. 
cd ('data/genes/processedData')
load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq%dDistThresh%d.mat', parcellation, distThreshold))
cd ..
cd 'parcellations'
[~, data_parcel]=read('MMPinMNI.nii');
ROIs = 1:180; 
coordsAssigned = samples2MNIparcellation(DataCoordinatesMNI, data_parcel, ROIs); 
shortestDist = distanceONsurfaceMNI(data_parcel, ROIs, coordsAssigned); 

