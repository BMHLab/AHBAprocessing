%------------------------------------------------------------------------------
% This script defines options to create ROI x gene matrices for each of the 4 brain parcellations
% First, create the initial data:
%------------------------------------------------------------------------------
% 1. import otiginal data from excel to matlab
% 2. perform probe selection

clear all;
close all;

options.ExcludeCBandBS = true;
options.useCUSTprobes = true; 
options.updateProbes = 'reannotator';
options.probeSelections = {'RNAseq'};
options.signalThreshold = 0.5;
options.RNAseqThreshold = 0.2;
options.RNAsignThreshold = false;
options.saveOutput = true;

S1_extractData(options)
S2_probes(options)

% 3. assign samples to parcellations, normalise expression measures and make the final data
%------------------------------------------------------------------------------
% Define general options that will be used in all cases
%------------------------------------------------------------------------------
options.correctDistance = false;
options.calculateDS = true;
options.distanceCorrection = 'Euclidean';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex';% what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.normaliseWithinSample = true;
options.xrange = [0 200];
options.plotCGE = false;
options.plotResiduals = false;
options.meanSamples = 'meanSamples'; %'meanSamples'; meanSubjects

%------------------------------------------------------------------------------
% 34 ROIs per hemisphere parcellation
%------------------------------------------------------------------------------

options.parcellations = {'aparcaseg'};
% for each parcellation first run with options.distanceThreshold = 40 to
% evaluate all distances between all assigned samples
options.distanceThreshold = 40;
S3_samples2parcellation(options)
% then use the appropriate threshold
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)

%------------------------------------------------------------------------------
% 100 ROIs per hemisphere parcellation
%------------------------------------------------------------------------------

options.parcellations = {'cust100'};
options.distanceThreshold = 40;
S3_samples2parcellation(options)
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)

%------------------------------------------------------------------------------
% 180 ROIs per hemisphere parcellation
%------------------------------------------------------------------------------
options.parcellations = {'HCP'};
options.distanceThreshold = 40;
S3_samples2parcellation(options)
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)

%------------------------------------------------------------------------------
% 250 ROIs per hemisphere parcellation
%------------------------------------------------------------------------------
options.parcellations = {'cust250'};
options.distanceThreshold = 40;
S3_samples2parcellation(options)
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)
