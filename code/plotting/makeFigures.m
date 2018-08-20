% % ------------------------------------------------------------------------------
% This script generates the relevant data before reproducing the figures.
% Some data that required long computation hours or different software are pre-calculated.
% At each step sequential data generation steps will be applied as needed
% % ------------------------------------------------------------------------------
% Generate initial data
options = struct();
options.ExcludeCBandBS = true;
options.useCUSTprobes = true;
options.updateProbes = 'reannotator';
fprintf('Generating the initial data\n')
S1_extractData(options)

% % ------------------------------------------------------------------------------
% Reproducing the figures
% % ------------------------------------------------------------------------------
% Figure 3
% % ------------------------------------------------------------------------------
figure3()

% % ------------------------------------------------------------------------------
% Figure 4AC
% % ------------------------------------------------------------------------------
options.saveOutput = true;
data2use = 'precomputed data';
% 'generate data': all data will be generated from original files - this will take >12h on a PC;
% 'precomputed data': correlation matrices will be loaded from pre-computed files.
afterQC = false;
% true: calculations will be performed using data after intensity-based
% filtering applying 0.5 threshold (Figure S2)
% false: calculations will be performed on original data without intensity
% based filtering ((Figure 4AC)
Regeneratedata = false;
fprintf('Making figure 4 A and C\n')
figure4AC(data2use, afterQC, options, Regeneratedata)

% % ------------------------------------------------------------------------------
% Figure 4B
% % ------------------------------------------------------------------------------
options.probeSelections = {'RNAseq'};
options.signalThreshold = -1; % this indicated no intensity based filtering
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
fprintf('Generating initial data for figure 4B\n')
S2_probes(options);
fprintf('Making figure 4B\n')
figure4B()
% % ------------------------------------------------------------------------------
% Figure 5C
% % ------------------------------------------------------------------------------
options.distanceThreshold = 40;
options.parcellations = {'aparcaseg'};
fprintf('Generating initial data for figure 5C\n')
S3_samples2parcellation(options)
options.distanceThreshold = 2;
S3_samples2parcellation(options)
fprintf('Making figure 5C\n')
figure5C()
% % ------------------------------------------------------------------------------
% Figure 6A
% % ------------------------------------------------------------------------------
% create data with intensity-based filtering
options.probeSelections = {'RNAseq'};
options.signalThreshold = 0.5; % this indicated no intensity based filtering
options.RNAseqThreshold = 0.2; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
fprintf('Generating initial data for figure 6\n')
S2_probes(options);

options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
S3_samples2parcellation(options)

doNormalise = false;
numROIs = 34;
% 34 - for left cortex only
% 41 - for left cortex + subcortex
fprintf('Making figure 6A\n')
figure; set(gcf,'Position',[300 300 1000 1000])
subplot(2,2,1); figure6ABCD(doNormalise,numROIs)
% % ------------------------------------------------------------------------------
% Figure 6B
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'zscore';
fprintf('Making figure 6B\n')
subplot(2,2,2); figure6ABCD(doNormalise,numROIs,whatNormalisation)

% % ------------------------------------------------------------------------------
% Figure 6C
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'scaledRobustSigmoid';
fprintf('Making figure 6C\n')
subplot(2,2,3); figure6ABCD(doNormalise,numROIs,whatNormalisation)

% % ------------------------------------------------------------------------------
% Figure 6D
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'scaledRobustSigmoid';
uselimma = true;
fprintf('Making figure 6D\n')
subplot(2,2,4); figure6ABCD(doNormalise,numROIs,whatNormalisation, uselimma)

% % ------------------------------------------------------------------------------
% Figure 6E
% % ------------------------------------------------------------------------------
numROIs = 34;
fprintf('Making figure 6E\n')
figure6E(numROIs)

% % ------------------------------------------------------------------------------
% Figure 8D
% % ------------------------------------------------------------------------------
fprintf('Making figure 8D\n')
figure8D()

% % ------------------------------------------------------------------------------
% Figure 9AB
% % ------------------------------------------------------------------------------
options.correctDistance = true;
options.calculateDS = true;
options.distanceCorrection = 'Surface';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = false;
options.xrange = [20 220];
options.plotCGE = true;
options.plotResiduals = true;
options.meanSamples = 'meanSamples';

fprintf('Generating data and making figure 9\n')
c = S4_normalisation(options)
% % ------------------------------------------------------------------------------
% Figure S1
% % ------------------------------------------------------------------------------
fprintf('Making figure S1\n')
figureS1()


% % ------------------------------------------------------------------------------
% Figure S2
% % ------------------------------------------------------------------------------
data2use = 'precomputed data';
afterQC = true;
Regeneratedata = false;
fprintf('Making figure S2\n')
figure4AC(data2use, afterQC, options, Regeneratedata)
pause(5)

% % ------------------------------------------------------------------------------
% Figure S3
% % ------------------------------------------------------------------------------
% A)
doNormalise = false;
numROIs = 41;
% 34 - for left cortex only
% 41 - for left cortex + subcortex
fprintf('Making figure S3\n')
figure; set(gcf,'Position',[300 300 1400 500])
subplot(1,2,1);
figure6ABCD(doNormalise,numROIs)

% B)
doNormalise = true;
whatNormalisation = 'scaledRobustSigmoid';
numROIs = 41;
% 34 - for left cortex only
% 41 - for left cortex + subcortex
subplot(1,2,2);
figure6ABCD(doNormalise,numROIs)

% % ------------------------------------------------------------------------------
% Figure S5
% % ------------------------------------------------------------------------------
% A,B)
options.correctDistance = true;
options.calculateDS = true;
options.distanceCorrection = 'GMvolume';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = false;
options.xrange = [20 190];
options.plotCGE = true;
options.plotResiduals = true;

fprintf('Making figure S4 A and B \n')
c = S4_normalisation(options)

% C)D)
options.correctDistance = true;
options.calculateDS = true;
options.distanceCorrection = 'Euclidean';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = false;
options.xrange = [20 160];
options.plotCGE = true;
options.plotResiduals = true;

fprintf('Making figure S4 C and D \n')
c = S4_normalisation(options)

% % ------------------------------------------------------------------------------
% Figure S6
% % ------------------------------------------------------------------------------
% generate relevant data files
options.correctDistance = false;
options.calculateDS = true;
options.distanceCorrection = 'SurfaceANDEuclidean';
options.Fit = {'exp'};
options.normaliseWhat = 'LcortexSubcortex'; %'LcortexSubcortexSEPARATE';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = false;
options.xrange = [0 220];
options.plotCGE = false;
options.plotResiduals = false;

fprintf('Generating data for figure S6A \n')
S4_normalisation(options);

howTOnormalise = 'together';
fprintf('Making figure S6A \n')
figureS6(howTOnormalise)

options.normaliseWhat = 'LcortexSubcortexSEPARATE'; %'LcortexSubcortexSEPARATE';% what part of the brain is normalised
fprintf('Generating data for figure S6B \n')
S4_normalisation(options);

howTOnormalise = 'separately';
fprintf('Making figure S6B \n')
figureS6(howTOnormalise)

% % ------------------------------------------------------------------------------
% Figure S7
% % ------------------------------------------------------------------------------
% generate relevant data files
options.correctDistance = true;
options.calculateDS = true;
options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
options.distanceCorrection = 'Euclidean';
options.Fit = {'removeMean'};
options.normaliseWhat = 'Lcortex'; %'LcortexSubcortexSEPARATE';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = false;
options.plotCGE = true;
options.plotResiduals = false;
options.percentDS =  10;
options.xrange = [20 160];
% A) 10% DS genes, low resolution parcellation
fprintf('Making figure S7A \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'low wesolution parcellation'});

% B) 100% DS genes, low resolution parcellation
options.percentDS =  100;
fprintf('Making figure S7B \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'low wesolution parcellation'});

% C) 10% DS genes, high resolution parcellation
options.parcellations = {'HCP'};
options.distanceThreshold = 40;
S3_samples2parcellation(options)
options.distanceThreshold = 2;
S3_samples2parcellation(options)
options.percentDS =  10;
options.xrange = [0 160];
fprintf('Making figure S7C \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'high wesolution parcellation'});

% D) 100% DS genes, high resolution parcellation
options.percentDS =  100;
options.xrange = [0 160];
fprintf('Making figure S7D \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'high wesolution parcellation'});
