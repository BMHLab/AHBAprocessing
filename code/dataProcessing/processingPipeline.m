
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Selects a single probe to represent a gene
%   3. assigns samples to the parcellations for each subject
%   4. normalises gene expression measures
%%
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
clear all;
close all;

options = struct();
options.ExcludeCBandBS = true;
options.useCUSTprobes = true;
options.updateProbes = 'reannotator';
options.probeSelections = {'RNAseq'};
options.parcellations = {'HCP'};
options.distanceThreshold = 2;
options.signalThreshold = 0.5;
options.RNAseqThreshold = 0.2;
options.RNAsignThreshold = false;
options.correctDistance = false;
options.calculateDS = true;
options.distanceCorrection = 'Euclidean';
options.Fit = {'exp'};
options.normaliseWhat = 'Lcortex'; % what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.resolution = 'ROI';
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [0 200];
options.plotCGE = true;
options.plotResiduals = true;
options.meanSamples = 'meanSamples'; %'meanSamples'; meanSubjects

%-------------------------------------------------------------------------------
% Extract data from excel files and saves to Matlab format:
S1_extractData(options)
%-------------------------------------------------------------------------------
S2_probes(options)
% for each parcellation first run with options.distanceThreshold = 40;
options.distanceThreshold = 40;
S3_samples2parcellation(options)
% then use the appropriate threshold
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)
