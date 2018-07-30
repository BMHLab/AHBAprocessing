This is a README file for the processingPipeline.m script

%------------------------------------------------------------------------------
General comments
%------------------------------------------------------------------------------
1. Probe-to-gene re-annotation provided with this data has been performed at the beginning of 2018. When performing analyses, it is recommended to use the latest available data, therefore new probe-to-gene re-annotations should be performed. This can be done using different software packages:
	- Re-annotator (software: https://sourceforge.net/projects/reannotator/; article: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139516). For more information on how to perform the re-annotation, please see README_Reannotator.txt
	- biomart (http://asia.ensembl.org/biomart/martview/64a0f53ea90ed52bf1641929d9fe7f9f)
	- AnnotationDbi (https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html).
2. Distances between samples on the cortical surface and within grey matter volume were calculated using aparcaseg parcellation (34 regions per hemisphere + 7 subcortical regions). Samples were matched to the regions of interest by generating parcellations in subject-specific spaces for every brain, therefore the number of mapped samples (after applying a selected distance threshold) onto each parcellation is slightly different. In order to use distances between samples calculated within the grey matter volume and on the cortical surface for other parcellations (HCP, cust100 or cust250), the mapping and calculations should be performed separately.
For more details on how this was done, please refer to:
	- calculateDistancesGMinMNI.m, a script that contains functions to calculate distances between samples within the grey matter volume.
	- dist_S1…, dist_S2…, dist_S3…, scripts to calculate distances on the cortical surface. This involves using freesurfer - the example of code is provided in runSurf2surf.sh script.

%------------------------------------------------------------------------------
Comments about processingPipeline.m script
%------------------------------------------------------------------------------
This script:
  1. Allows to specify options for the data processing including:
        - importing data from original excel sheets to matlab files: S1_extractData.m
        - selecting a representative probe for a gene: S2_probes.m
        - assigning samples to parcellations: S3_samples2parcellation.m
        - normalizing the data, making region x gene matrix and CGE matrix: S4_normalisation.m

%------------------------------------------------------------------------------
To choose the options specify selections in the following manner.
%------------------------------------------------------------------------------
options.ExcludeCBandBS = true;
true: Brainstem and cerebellum samples will be excluded based on labels provided by the AHBA
false: Brainstem and cerebellum samples will be retained.
%------------------------------------------------------------------------------
options.useCUSTprobes = true;
true: CUST and Agilent probes will be used
false: only Agilent probes will be used
%------------------------------------------------------------------------------
options.updateProbes = 'reannotator';
this pipeline does not perform the re-annotation, it should be done
separately and corresponding files need to be used. If new files are
being made, functions to read them into matlab should be updated/checked.

'Biomart': probe-to-gene annotation file is generated using Biomart web-based tool:
http://asia.ensembl.org/biomart/martview/e4be7e990dff3b5ebd5070366bc81c8c;
filename: mart_export_updatedProbes.txt located in data/genes/rawData

'reannotator': probe-to-gene annotation file is generated using Re-annotator:
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139516
filename: reannotatedProbes.mat located in data/genes/rawData

'no': no probe-to-gene annotation is performed. Data presented by AHBA is
used.
%------------------------------------------------------------------------------
options.probeSelections = {'RNAseq'};
probe selection is performed for each subject separately. Then the probe
that is chosen in most subjects is selected to represent the gene for all
subjects.

'Variance': probe with the highest variance across regions is selected

'RNAseq': probe with the highest correlation to RNAseq data is selected

'PC': probe with the highest loading on the 1st principal component is
selected. No z-scoring is performed on the data to reproduce the result
in Parkes, 2017.

'maxIntensity': probe with the highest mean expression value is selected.

'maxCorrelation_intensity': probe with the highest correlation to other
probes is selected when >2 probes are available. When 2 probes are
available, probe with the higher average intensity is selected.

'maxCorrelation_variance': probe with the highest correlation to other
probes is selected when >2 probes are available. When 2 probes are
available, probe with the higher variance is selected.

'CV': probe with the highest coefficient of variation across regions is selected

'LessNoise': probe with the highest proportion of samples with expression values exceeding the background is selected.
Binary indication of expression values exceeding the background are
determined by AHBA.

'Mean': mean of all availabel probes for a gene is calculated.

'Random': a probe is selected at random.

'DS': probe with the higher differential stability score is selected.
%------------------------------------------------------------------------------
options.parcellations = {'HCP'};
'HCP': 180 nodes per hemisphere
'aparcaseg': 34 nodes per hemisphere + subcortex
'cust100': 100 nodes per hemisphere + subcortex
'cust250': 250 nodes per hemisphere + subcortex
%------------------------------------------------------------------------------
options.distanceThreshold = 2;
distance threshold imposed when assigning samples to regions of the parcellation
%------------------------------------------------------------------------------
options.signalThreshold = 0.5;
the proportion of samples with expression values exceeding the
background.
%------------------------------------------------------------------------------
options.RNAseqThreshold = 0.2;
the minimal correlation between microarray and RNA-seq expression
measures required. Genes with correlation values lower then this
threshold will be excluded.
%------------------------------------------------------------------------------
options.RNAsignThreshold = false;
false: to significance threshold on RNAseq vs microarray correlation is required;
true: 0.05 significance threshold for RNAseq vs microarray correlation is required;
%------------------------------------------------------------------------------
options.correctDistance = false;
false: no CGE-distance correction is performed
true: CGE-distance correction is performed, refer to
options.distanceCorrection for more options;
%------------------------------------------------------------------------------
options.calculateDS = true;
true: evaluate the differential stability for each gene
false: do not evaluate the differential stability for each gene
%------------------------------------------------------------------------------
options.distanceCorrection = 'Euclidean';
'Surface': distances between samples are evaluated on the
corticalsurface. Values are pre-computed for aparcaseg (34 nodes per
hemisphere) parcellation. Calculating distances on surface required
additional software.

'GMvolume': distances between samples are evaluated within the cortical
grey matter volume. Values are pre-computed for aparcaseg (34 nodes per
hemisphere) parcellation. Calculation requires time.

'Euclidean': distances between samples are evaluated as straight lines
using Euclidean distances between sample coordinates. Can be performed on
all parcellations including cortical and subcortical regions.

'SurfaceANDEuclidean': distances between cortical samples are evaluated
on the cortical surface for aparcaseg (34 nodes perhemisphere) parcellation.
Distances between subcortical regions as well as between cortical and
subcortical regions evaluated as straight lines using Euclidean distance.
%------------------------------------------------------------------------------
options.Fit = {'exp'};
'exp': exponential fit to approximate CGE vs distance
relationship:  FitCurve = c.A*exp(-c.n*Dvect) + c.B;
'removeMean': the mean of each bin is removed to approximate CGE vs distance relationship
'linear': linear fit to approximate CGE vs distance relationship: FitCurve = c.p1*Dvect + c.p2;
'exp_1_0': 1-0 constrained exponential fit to approximate CGE vs distance
relationship: FitCurve = exp(-c.n*Dvect);
'exp0': 0 constrained exponential fit to approximate CGE vs distance
relationship: FitCurve = c.A.*exp(-c.n*Dvect);
'exp1': 1 constrained exponential fit to approximate CGE vs distance
relationship: FitCurve = exp(-c.n*Dvect) + c.B;
'decay': FitCurve = c.A/Dvect + c.B;
%------------------------------------------------------------------------------
options.normaliseWhat = 'Lcortex';what part of the brain is normalised
'Lcortex': only data from left cortex are used (1:6 subjects)
'LcortexSubcortex': only data from left cortex and subcortex are used,
normalised together (1:6 subjects)
'wholeBrain': data from the whole brain are used (cortex+ subcortex both
hemispheres), normalised together (1:2 subjects) given only 2 subjects
have data from the right hemisphere;
'LRcortex': data from the whole brain are used (cortex both
hemispheres), normalised together (1:2 subjects) given only 2 subjects
have data from the right hemisphere;
'LcortexSubcortexSEPARATE': only data from left cortex and subcortex are used,
normalised separately (1:6 subjects)
%------------------------------------------------------------------------------
options.normMethod = 'scaledRobustSigmoid'; what type of normalisation method used
'subtractMean'Subtract the mean:
'maxmin' Linear rescaling to the unit interval
'zscore'
'robustSigmoid'A outlier-robust sigmoid
'scaledRobustSigmoid' A scaled, outlier-robust sigmoid
'sigmoid' Standard sigmoidal transformation
'scaledSigmoid' Standard sigmoid transform, then a rescaling to the unit interval
'mixedSigmoid' Uses a scaled sigmoid if iqr=0; a scaled, outlier-robust sigmoid otherwise; Uses only non-NaN; Outlier-adjusted sigmoid:
'scaledsigmoid5q' First caps at 5th and 95th quantile, then does scaled sigmoid
%------------------------------------------------------------------------------
options.percentDS =  100;
the percentage of highest differential stability genes used;
%------------------------------------------------------------------------------
options.doNormalise = true;
true: normalise gene expression (each subject separately)
false: do not normalise gene expression
%------------------------------------------------------------------------------
options.resolution = 'ROI';
'ROI': expression measures are first averaged into ROIs and then CGE is
calculated
'sample': CGE is calculated at the resolution of samples and then CGE
values are averaged to corresponding pairs of ROIs
%------------------------------------------------------------------------------
options.saveOutput = true;
%------------------------------------------------------------------------------
options.normaliseWithinSample = false;
true: normalise expression values within sample before within gene normalisation
false: do not normalise expression values within sample before within
gene normalisation
%------------------------------------------------------------------------------
options.meanSamples = 'meanSamples';
'meanSamples' when multiple samples are available in a region, the mean of all samples is calculated to summarise the expression vector. The “weight” of every sample in that region is equal.
‘meanSubjects’ when multiple samples are available in a region, the average of samples for each subject is calculated first, and then the mean expression determined as an average across the brains. The “weight” of every subjects that has samples in that region is equal.
see Figure 7 for more details.
%------------------------------------------------------------------------------

S1_extractData(options)
S2_probes(options)
for each parcellation first run with options.distanceThreshold = 40;
options.distanceThreshold = 40;
S3_samples2parcellation(options)
then use the appropriate threshold
options.distanceThreshold = 2;
S3_samples2parcellation(options)
c = S4_normalisation(options)
