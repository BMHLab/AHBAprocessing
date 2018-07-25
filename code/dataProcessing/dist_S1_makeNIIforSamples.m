
%% This script:
%   1. Makes a .nii file for all samples in a brain for each subject
%%
%------------------------------------------------------------------------------
% for each subject make a .nii file with all zeros and assign indexes for the coordinates of mappped samples.
%------------------------------------------------------------------------------

subjects = 1:6;
brainPart = 'Lcortex';
numSamples = 1;
parcellation = {'aparcaseg'};%, 'cust100', 'cust250'};
distanceThreshold = 2;
sampleIND = cell(6,1);
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqQC82DistThresh2.mat')
cd ../../..
for subject = subjects
    cd ('data/genes/parcellations')
    subjectDir = sprintf('S0%d_H0351', subject);
    cd (subjectDir)

    fprintf('Subject %u parcellation %s assignment distance threshold %u\n; ', subject, parcellation{1}, distanceThreshold )

    %------------------------------------------------------------------------------
    % Load parcellations
    %------------------------------------------------------------------------------
    if strcmp(parcellation, 'aparcaseg')
        parcName = 'default_NativeAnat';
        %cd (parcName);
        [hdr, data_parcel]=read('defaultparc_NativeAnat.nii');
        NumNodes = 82;
        LeftCortex = 1:34;
        LeftSubcortex = 35:41;
        RightCortex = 42:75;
        RightSubcortex = 76:82;
    elseif strcmp(parcellation, 'cust100')
        parcName = 'custom100_NativeAnat';
        %cd (parcName);
        [hdr, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 220;
        LeftCortex = 1:100;
        LeftSubcortex = 101:110;
        RightCortex = 111:210;
        RightSubcortex = 211:220;
    elseif strcmp(parcellation, 'cust250')
        parcName = 'custom250_NativeAnat';
        %cd (parcName);
        [hdr, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 530;
        LeftCortex = 1:250;
        LeftSubcortex = 251:265;
        RightCortex = 266:515;
        RightSubcortex = 516:530;

    end
    % make an empty image
    image = zeros(size(data_parcel));
    coordinates = DataCoordinatesMRI{subject}(:,3:5);
    switch brainPart
        case 'Lcortex'
            nROIs = LeftCortex;
        case 'LcortexSubcortex'
            nROIs = 1:max(LeftSubcortex);
        case 'wholeBrain'
            nROIs = 1:NumNodes;
        case 'LRcortex'
            nROIs = [LeftCortex,RightCortex];
    end
    samplesIND = find(ismember(DataCoordinatesMRI{subject}(:,2),nROIs));
    coordinates2 = coordinates(samplesIND,:);
    index = numSamples:numSamples+max(samplesIND)-1;

    [C,ia,ic] = unique(coordinates2, 'rows', 'stable');
    numSamples = numSamples+length(ia);
    sampleIND{subject} = ia;
    coordinates3 = coordinates2(ia,:);
    for samp=1:length(ia)
        point = coordinates3(samp,:);
        image(point(1), point(2), point(3)) = index(samp);
    end
    cd ../../..
    cd ('genes/processedData')
    filename1 = sprintf('S%dsamplesXXX.nii', subject);
    write(hdr,image,filename1);

    cd ../../..

end
cd ('data/genes/processedData')
% when calculating sample-smple coexpression keep only those samples
save('samples2keepXXX.mat', 'sampleIND');

%% after running this script
% 1. we need to transform volume files (.nii) to surface files (.mgz) using mri_vol2surf command on massive (M2: /gpfs/M2Scratch/Monash076/aurina/AllenInstitute) for each subject and produce one .mgz file per subject
% mri_vol2surf --mov S1samples.nii --regheader H0351_2001 --hemi lh --o S1samplesTEST.mgz --projfrac-max 0 1 0.1 --interp nearest
