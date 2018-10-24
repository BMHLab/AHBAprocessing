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
%S1_extractData(options)

% % ------------------------------------------------------------------------------
% Reproducing the figures
% % ------------------------------------------------------------------------------
% % ------------------------------------------------------------------------------
% Figure 1 - data schematic
% % ------------------------------------------------------------------------------
% % ------------------------------------------------------------------------------
% Figure 2 - processing schematic
% % ------------------------------------------------------------------------------
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
% filtering applying 0.5 threshold (Figure S1)
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
options.VARfilter = false;
options.VARscale = 'log2';
options.VARperc = 50;
QClabel = 'noQC';
fprintf('Generating initial data for figure 4B\n')
S2_probes(options);
fprintf('Making figure 4B\n')
figure4B(QClabel)
% % ------------------------------------------------------------------------------
% Figure 5C
% % ------------------------------------------------------------------------------
options.probeSelections = {'RNAseq'};
options.signalThreshold = -1; % this indicated no intensity based filtering
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
options.VARfilter = false;
options.divideSamples = 'listCortex';
options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
options.excludeHippocampus = false;
fprintf('Generating initial data for figure 5C\n')
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
options.VARfilter = false;
options.divideSamples = 'listCortex';
options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
options.excludeHippocampus = false;
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
subplot(2,2,1); figure6ABCD(doNormalise,numROIs);
% % ------------------------------------------------------------------------------
% Figure 6B
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'zscore';
fprintf('Making figure 6B\n')
subplot(2,2,2); figure6ABCD(doNormalise,numROIs,whatNormalisation);

% % ------------------------------------------------------------------------------
% Figure 6C
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'scaledRobustSigmoid';
fprintf('Making figure 6C\n')
subplot(2,2,3); figure6ABCD(doNormalise,numROIs,whatNormalisation);

% % ------------------------------------------------------------------------------
% Figure 6D
% % ------------------------------------------------------------------------------
doNormalise = true;
numROIs = 34;
whatNormalisation = 'scaledRobustSigmoid';
uselimma = true;
fprintf('Making figure 6D\n')
subplot(2,2,4); figure6ABCD(doNormalise,numROIs,whatNormalisation, uselimma);

% % ------------------------------------------------------------------------------
% Figure 6E
% % ------------------------------------------------------------------------------
numROIs = 34;
fprintf('Making figure 6E\n')
figure6E(numROIs)

% % ------------------------------------------------------------------------------
% Figure 7D
% % ------------------------------------------------------------------------------
fprintf('Making figure 7D\n')
figure7D()

% % ------------------------------------------------------------------------------
% Figure 8AB
% % ------------------------------------------------------------------------------
options.correctDistance = true;
options.calculateDS = true;
options.distanceCorrection = 'Surface';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [20 220];
options.plotCGE = true;
options.plotResiduals = true;
options.meanSamples = 'meanSamples';
options.probeSelections = {'RNAseq'};
options.signalThreshold = 0.5; % this indicated no intensity based filtering
options.RNAseqThreshold = 0.2; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
options.VARfilter = false;
options.divideSamples = 'listCortex';
options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
options.excludeHippocampus = false;
fprintf('Generating initial data for figure 8\n')
S2_probes(options);
S3_samples2parcellation(options);

fprintf('Generating data and making figure 8\n')
c = S4_normalisation(options)
% % ------------------------------------------------------------------------------
% Figure S1
% % ------------------------------------------------------------------------------

data2use = 'precomputed data';
afterQC = true;
Regeneratedata = false;
fprintf('Making figure S1\n')
figure4AC(data2use, afterQC, options, Regeneratedata)
pause(5)

options.probeSelections = {'RNAseq'};
options.signalThreshold = 0.5; % this indicated intensity based filtering
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
QClabel = 'QC';
fprintf('Generating initial data for figure S2B\n')
S2_probes(options);
fprintf('Making figure S1B\n')
figure4B(QClabel)
% % ------------------------------------------------------------------------------
% Figure S2
% % ------------------------------------------------------------------------------
options.correctDistance = false;
options.calculateDS = true;
options.distanceCorrection = 'Euclidean';
options.Fit = {'exp'};
options.signalThreshold = -1; % this indicated no intensity based filtering
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.RNAsignThreshold = false;
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [20 190];
options.plotCGE = false;
options.plotResiduals = false;
S4_normalisation(options);
fprintf('Making figure S2\n')
[rDS, rRNA] = figureS2()

% % ------------------------------------------------------------------------------
% Figure S3
% % ------------------------------------------------------------------------------
% to compare to variance based filtering, generate data
% without the IBF and with VBF on normal-scaled data
options.signalThreshold = -1;
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.VARfilter = true;
options.VARscale = 'normal'; % 'log2', 'normal'
options.VARperc = 50; %
QClabel = sprintf('noQCvar%s', options.VARscale);
S2_probes(options);
figure4B(QClabel)
fprintf('Making figure S2\n')
% without the IBF and with VBF on log2-scaled data
options.signalThreshold = -1;
options.RNAseqThreshold = -1; % this indicated no exclusion of data based on correlation to RNA-seq
options.VARfilter = true;
options.VARscale = 'log2';
options.VARperc = 50; %
QClabel = sprintf('noQCvar%s', options.VARscale);
S2_probes(options);
figure4B(QClabel)

% % ------------------------------------------------------------------------------
% Figure S4
% % ------------------------------------------------------------------------------
fprintf('Making figure S4\n')
figureS4()
% % ------------------------------------------------------------------------------
% Figure S5
% % ------------------------------------------------------------------------------
% A)
doNormalise = false;
numROIs = 41;
% 34 - for left cortex only
% 41 - for left cortex + subcortex
fprintf('Making figure S5\n')
figure; set(gcf,'Position',[300 300 1400 500]);
subplot(1,2,1);
figure6ABCD(doNormalise,numROIs);

% B)
doNormalise = true;
whatNormalisation = 'scaledRobustSigmoid';
numROIs = 41;
% 34 - for left cortex only
% 41 - for left cortex + subcortex
subplot(1,2,2);
figure6ABCD(doNormalise,numROIs);

% % ------------------------------------------------------------------------------
% Figure S6 - normalization schematic
% % ------------------------------------------------------------------------------
% % ------------------------------------------------------------------------------
% Figure S7 - sample averaging schematic
% % ------------------------------------------------------------------------------
% % ------------------------------------------------------------------------------
% Figure S8
% % ------------------------------------------------------------------------------
% A,B)
options.signalThreshold = 0.5; % this indicated intensity based filtering
options.RNAseqThreshold = 0.2; % this indicated no exclusion of data based on correlation to RNA-seq
options.correctDistance = true;
options.calculateDS = true;
options.distanceCorrection = 'GMvolume';
options.Fit = {'exp'};
options.VARfilter = false;
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [0 200];
options.plotCGE = true;
options.plotResiduals = true;
S2_probes(options)
S3_samples2parcellation(options)
fprintf('Making figure S8 A and B \n')
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
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [0 200];
options.plotCGE = true;
options.plotResiduals = true;

fprintf('Making figure S8 C and D \n')
c = S4_normalisation(options)

% % ------------------------------------------------------------------------------
% Figure S9
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
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [0 220];
options.plotCGE = false;
options.plotResiduals = false;

fprintf('Generating data for figure S9A \n')
S4_normalisation(options);

howTOnormalise = 'together';
fprintf('Making figure S9A \n')
figureS9(howTOnormalise)

options.normaliseWhat = 'LcortexSubcortexSEPARATE'; %'LcortexSubcortexSEPARATE';% what part of the brain is normalised
fprintf('Generating data for figure S9B \n')
S4_normalisation(options);

howTOnormalise = 'separately';
fprintf('Making figure S9B \n')
figureS9(howTOnormalise)

% % ------------------------------------------------------------------------------
% Figure S10
% % ------------------------------------------------------------------------------
% generate relevant data files
options.correctDistance = true;
options.calculateDS = true;
options.parcellations = {'aparcaseg'};
options.distanceThreshold = 2;
options.distanceCorrection = 'Surface';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex'; %'LcortexSubcortexSEPARATE';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.doNormalise = true;
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.plotCGE = true;
options.plotResiduals = false;
options.percentDS =  10;
options.xrange = [0 220];
% A) 10% DS genes, low resolution parcellation
fprintf('Making figure S10A \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'low wesolution parcellation'});

% B) 100% DS genes, low resolution parcellation
options.percentDS =  100;
fprintf('Making figure S10B \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'low wesolution parcellation'});

% C) 10% DS genes, high resolution parcellation
options.parcellations = {'HCP'};
options.distanceThreshold = 2;
S3_samples2parcellation(options)
options.percentDS =  10;
options.xrange = [0 220];
fprintf('Making figure S10C \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'high wesolution parcellation'});

% D) 100% DS genes, high resolution parcellation
options.percentDS =  100;
options.xrange = [0 220];
fprintf('Making figure S10D \n')
S4_normalisation(options);
title({sprintf('%d%% DS genes', options.percentDS); 'high wesolution parcellation'});
