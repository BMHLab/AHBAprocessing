
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
options.parcellations = {'cust100'};
options.divideSamples = 'ontology'; %'listCortex'; % 'ontology' - GM volume distances are pre-calculated only for listCortex option; 
options.excludeHippocampus = false;
options.distanceThreshold = 2;
options.signalThreshold = 0.5; % -1 implies no IBF filtering
options.VARfilter = false; 
options.VARscale = 'normal'; % 'log2', 'normal'
options.VARperc = 50;
options.RNAsignThreshold = false;
options.RNAseqThreshold = 0.2;
options.correctDistance = false;
options.calculateDS = true;
options.distanceCorrection = 'Euclidean';
options.Fit = {'exp'};
options.normaliseWhat = 'LRcortex'; % what part of the brain is normalised
options.normMethod = 'scaledRobustSigmoid'; % what type of normalisation method used
options.percentDS =  100;
options.doNormalise = true;
options.saveOutput = true;
options.normaliseWithinSample = true;
options.xrange = [0 220];
options.plotCGE = true;
options.plotResiduals = true;
options.meanSamples = 'meanSamples'; %'meanSamples'; meanSubjects

%-------------------------------------------------------------------------------
% Extract data from excel files and saves to Matlab format:
S1_extractData(options)
%-------------------------------------------------------------------------------
S2_probes(options)
S3_samples2parcellation(options)
c = S4_normalisation(options)
