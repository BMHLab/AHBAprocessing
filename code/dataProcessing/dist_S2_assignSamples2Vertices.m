%% Author: Aurina
%% Date modified: 2017-08-18
%% This script:
%   1. loads .mgz files for each brain (with all samples in them)
%   2. loads original orig.mgz file, pial surface of the brain, original
%   001.nii file, and .mgz file that we created in the previous step using
%   mri_vol2surf for all samples
%   3. Transforms coodrinates
%   4. Creates a separate .mgz file per each sample and a list of samples
%   for each subject
%%

%------------------------------------------------------------------------------
% Load data where each sample is assigned to parcellation, so we know which
% sample belongs where. Aonly samples that were already assigned to the
% parcellation with 2mm threshold are chosen here
%------------------------------------------------------------------------------

cd ('data/genes/processedData')
%load('MicroarrayDatadPC82DistThresh2_CoordsAssigned.mat');
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqQC82DistThresh2.mat')
cd ../forFreesurfer
keepSamplesOrig = cell(6,1);
keepSamples = cell(6,1);

%------------------------------------------------------------------------------
% Load fsaverage surface template for the left pial surface
%------------------------------------------------------------------------------
[verticesFSaverage,facesFSaverage] = read_surf('lhfsaverage.white');

%------------------------------------------------------------------------------
% For each subject and each sample create a separate .mgz file
%------------------------------------------------------------------------------
k=1;
for sub = 1:6
    cd (sprintf('S0%d', sub))
    % Here first load up the freesurfer volume and lh pial surface
    data = MRIread('orig.mgz');
    [vertices,faces] = read_surf('lh.pial');
    
    vert2 = [vertices-1 ones(size(vertices,1),1)];
    h2 = inv(data.tkrvox2ras)*vert2.';
    
    data2 = MRIread('001.nii');
    h3 = inv(data.tkrvox2ras*inv(data.vox2ras)*(data2.niftihdr.sform))*vert2.';
    
    % Vertices in voxel space for samples and T1
    newVertices = h3(1:3,:).';
    [xx, yy, zz] = meshgrid(1:size(data2.vol,2),1:size(data2.vol,1),1:size(data2.vol,3));
    
    % load overlay for a subject
    dataOrig = MRIread(sprintf('S%dsamplesXXX.mgz', sub));
    
    coordinatesall = DataCoordinatesMRI{sub};
    ROI = coordinatesall(:,2); keep = find(ROI<=34);
    coordinates = coordinatesall(keep,3:5);
    
    coordinatesNEWvox = zeros(length(keep),3);
    coordinatesNEWvert = zeros(length(keep),3);
    overlay = zeros(size(vertices));
    %int=linspace(1,1249,1249);
    
    for i=1:length(keep)
        dataOrig.vol = zeros(size(dataOrig.vol));
    
        x = coordinates(i,1);
        y = coordinates(i,2);
        z = coordinates(i,3);
        
        [minval, ind] =min(sqrt((newVertices(:,1)-x).^2 + (newVertices(:,2)-y).^2 + (newVertices(:,3)-z).^2));
        vertexind = ind;
        
        % save coordinates
        coordinatesNEWvox(i,:) = [newVertices(vertexind,1),newVertices(vertexind,2),newVertices(vertexind,3)];
        coordinatesNEWvert(i,:) = [vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)];
        
        lab = i+k-1; 
        overlay(vertexind) = i;
        dataOrig.vol(vertexind) = i;
        
        MRIwrite(dataOrig,sprintf('S%dsample%d_singleVertXXX.mgz', sub, i));
        % In freesurfecoordinatesNEWvertr space vertex co-ordinate is:
        % disp([vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)]);
    end
    
    [~,ia] = unique(coordinatesNEWvert, 'rows', 'stable');
    keepSamplesOrig{sub} = keep;
    k=k+length(ia);
    keepSamples{sub} = ia;
    sampleList = keepSamplesOrig{sub}; 
    
%------------------------------------------------------------------------------
% For each subject create a list of subjects
%------------------------------------------------------------------------------
    fileID = fopen(sprintf('S%dsampleListXXX.txt', sub),'w');
    nbytes = fprintf(fileID,'%1d\n',sampleList); 
    fclose(fileID);
    cd ..
    
end
save('keepSamplesXXX.mat', 'keepSamples');
%% after running this script 
% 1. the folder forFreesurfer need to be copied to massive (M2: /gpfs/M2Scratch/Monash076/aurina/HumanExpression/data/genes/forFreesurfer)
% 2. then run runSurf2surf.sh script to map samples from native space on the surface to fsaverage standard space. 
% 3. This way for each sample a new .mgz file is created with multiple verteces assigned to one sample
% 4. Copy forFreesurfer folder back with all the necessary files 
% 5. Next script (dist_S3_chooseVertex4sample.m) will choose one vertex per sample and will calculate distances on the fsaverage surface
