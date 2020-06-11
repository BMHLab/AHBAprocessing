%% This script uses:
%   1. parcellation annotation in the fsaverage space (separate file should be used for different parcellations)
%   2. fsaverage file of the sphere
%   3. fsaverage file of the pial (or white matter) surface
%%
function distRegions = distanceONsurfaceROIs(parcellation, whatSurface)
% inputs:
% 1. parcellation: HCP, aparcaseg, cust100, cust250
% 2. whatSurface: pial, white
% Output:
% distances between regions on the surface for a selected parcellation
% load parcellation annotation file on the surface
if strcmp(parcellation, 'aparcaseg')
    [vertices, label,ctab] = read_annotation('lh.aparc.annot');
    LeftCortex = [2:4,6:36];
elseif strcmp(parcellation, 'cust100')
    [vertices, label,ctab] = read_annotation('lh.random200.annot');
    LeftCortex = 2:101;
elseif strcmp(parcellation, 'cust250')
    [vertices, label,ctab] = read_annotation('lh.random500.annot');
    LeftCortex = 2:251;
elseif strcmp(parcellation, 'HCP')
    [vertices, label,ctab] = read_annotation('lh.HCP-MMP1.annot');
    LeftCortex = 2:181;
elseif strcmp(parcellation, 'Schaefer100')
    [vertices, label,ctab] = read_annotation('lh.Schaefer100_7net.annot');
    LeftCortex = 2:51;
elseif strcmp(parcellation, 'Schaefer300')
    [vertices, label,ctab] = read_annotation('lh.Schaefer300_7net.annot');
    LeftCortex = 2:151;
elseif strcmp(parcellation, 'Schaefer500')
    [vertices, label,ctab] = read_annotation('lh.Schaefer500_7net.annot');
    LeftCortex = 2:251;
elseif strcmp(parcellation, 'Schaefer1000')
    [vertices, label,ctab] = read_annotation('lh.Schaefer1000_7net.annot');
    LeftCortex = 2:501;
end

% load sphere file
[verticesSphere] = read_surf('lh.sphere');
% load fsaverage surface
[verticesFSaverage,facesFSaverage] = read_surf(sprintf('lhfsaverage.%s', whatSurface));

%------------------------------------------------------------------------------
% For every region find vertices on the sphere and take the mean of their coordinates
%------------------------------------------------------------------------------

vertSurf = zeros(length(LeftCortex),1);
j=1;
for r=LeftCortex
    verts = find(label==ctab.table(r,5));
    coordROI = mean(verticesSphere(verts,:),1);
    % find the shortest distance between this center point and any vertex
    % from the region
    k = dsearchn(verticesSphere(verts,:),coordROI);
    % select that vertex in the parcellation
    vertSurf(j) = verts(k);
    j=j+1;
end

%------------------------------------------------------------------------------
% Calculate distances on the surface for each ROI
%------------------------------------------------------------------------------
distances = zeros(size(verticesFSaverage,1),length(vertSurf));

tic
for region=1:length(vertSurf)
    % calculate the distance between the selected vertex and all the
    % vertices in the parcellation
    [D] = perform_fast_marching_mesh(verticesFSaverage.', facesFSaverage.' + 1, vertSurf(region));
    distances(:,region) = D;
    fprintf('%d Region has been analised\n', region)
end
% select distances between the vertices that were representative of the ROI
distRegions = distances(vertSurf,:);
end
