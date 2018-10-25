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

'Mean': mean of all available probes for a gene is calculated.

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
options.divideSamples = 'listCortex';
'ontology' - samples are separated into cortical and subcortical based on AHBA ontology information
'listCortex' - samples are separated into cortical and subcortical based on re-defined list of region names;
%------------------------------------------------------------------------------
options.excludeHippocampus = false;
false - hippocampal samples are not excluded from sample assignment;
true - hippocampal samples are excluded from sample assignment;
%------------------------------------------------------------------------------
These options correspond to variance-based filtering
options.VARfilter = false;
false - variance based filtering is not performed; if this option is selected, other options such as options.VARscale and options.VARperc are not considered.
true - variance based filtering is performed;
options.VARscale = 'normal';
'normal' - variance based filtering performed on non-log2 transformed data (gives better outcomes)
'log2' - variance based filtering performed on log2 transformed data (gives worse outcomes)
options.VARperc = 50;
The percentage of lowest variance probes considered for filtering.
The probe is excluded if it is classified as low variance in all 6 subjects.
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
'Euclidean': distances between samples are evaluated as straight lines
using Euclidean distances between sample coordinates. Can be performed on
all parcellations including cortical and subcortical regions.

'Surface': distances between samples are evaluated on the cortical surface.
In a current configuration can be performed when evaluating expression only within left cortex. 

'GMvolume': distances between samples are evaluated within the cortical
grey matter volume. Values are pre-computed for aparcaseg (34 nodes per
hemisphere) parcellation. Calculation requires time.

'SurfaceANDEuclidean': distances between cortical regions are evaluated
on the cortical surface. Distances between subcortical regions as well as between cortical and
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
options.normaliseWhat = 'Lcortex'; what part of the brain is normalised
'Lcortex': only data from left cortex are used (1:6 subjects)
'LcortexSubcortex': only data from left cortex and subcortex are used,
normalised together (1:6 subjects)
'wholeBrain': data from the whole brain are used (cortex+ subcortex both
hemispheres), normalised together (1:2 subjects) given only 2 subjects
have data from the right hemisphere; this options is not suitable for HCP parcellation as it doesn't include subcortical regions.
'LRcortex': data from the whole brain are used (cortex both
hemispheres), normalised together (1:2 subjects) given only 2 subjects
have data from the right hemisphere;
'LcortexSubcortexSEPARATE': only data from left cortex and subcortex are used,
normalised separately (1:6 subjects)
%------------------------------------------------------------------------------
options.normMethod = 'scaledRobustSigmoid'; what type of normalisation method used
'subtractMean'Subtract the mean:
'maxmin' - Linear rescaling to the unit interval
'zscore' - z-score normalisation;
'robustSigmoid' - outlier-robust sigmoid
'scaledRobustSigmoid' - scaled, outlier-robust sigmoid
'sigmoid' - standard sigmoidal transformation
'scaledSigmoid' - standard sigmoid transform, then a rescaling to the unit interval
'mixedSigmoid' - uses a scaled sigmoid if iqr=0; a scaled, outlier-robust sigmoid otherwise; Uses only non-NaN; Outlier-adjusted sigmoid:
'scaledsigmoid5q' - first caps at 5th and 95th quantile, then does scaled sigmoid
%------------------------------------------------------------------------------
options.percentDS =  100;
the percentage of highest differential stability genes used;
%------------------------------------------------------------------------------
options.doNormalise = true;
true: normalise gene expression (each subject separately)
false: do not normalise gene expression
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
see Figure S7 for more details.
%------------------------------------------------------------------------------

S1_extractData(options)
S2_probes(options)
S3_samples2parcellation(options)
c = S4_normalisation(options)
%------------------------------------------------------------------------------
Additional information
%------------------------------------------------------------------------------
Region labels for aparcaseg parcellation
1    ctx-lh-bankssts
2    ctx-lh-caudalanteriorcingulate
3    ctx-lh-caudalmiddlefrontal
4    ctx-lh-cuneus
5    ctx-lh-entorhinal
6    ctx-lh-fusiform
7    ctx-lh-inferiorparietal
8    ctx-lh-inferiortemporal
9    ctx-lh-isthmuscingulate
10    ctx-lh-lateraloccipital
11    ctx-lh-lateralorbitofrontal
12    ctx-lh-lingual
13    ctx-lh-medialorbitofrontal
14    ctx-lh-middletemporal
15    ctx-lh-parahippocampal
16    ctx-lh-paracentral
17    ctx-lh-parsopercularis
18    ctx-lh-parsorbitalis
19    ctx-lh-parstriangularis
20    ctx-lh-pericalcarine
21    ctx-lh-postcentral
22    ctx-lh-posteriorcingulate
23    ctx-lh-precentral
24    ctx-lh-precuneus
25    ctx-lh-rostralanteriorcingulate
26    ctx-lh-rostralmiddlefrontal
27    ctx-lh-superiorfrontal
28    ctx-lh-superiorparietal
29    ctx-lh-superiortemporal
30    ctx-lh-supramarginal
31    ctx-lh-frontalpole
32    ctx-lh-temporalpole
33    ctx-lh-transversetemporal
34    ctx-lh-insula
35    Left-Thalamus-Proper
36    Left-Caudate
37    Left-Putamen
38    Left-Pallidum
39    Left-Hippocampus
40    Left-Amygdala
41    Left-Accumbens-area
42    ctx-rh-bankssts
43    ctx-rh-caudalanteriorcingulate
44    ctx-rh-caudalmiddlefrontal
45    ctx-rh-cuneus
46    ctx-rh-entorhinal
47    ctx-rh-fusiform
48    ctx-rh-inferiorparietal
49    ctx-rh-inferiortemporal
50    ctx-rh-isthmuscingulate
51    ctx-rh-lateraloccipital
52    ctx-rh-lateralorbitofrontal
53    ctx-rh-lingual
54    ctx-rh-medialorbitofrontal
55    ctx-rh-middletemporal
56    ctx-rh-parahippocampal
57    ctx-rh-paracentral
58    ctx-rh-parsopercularis
59    ctx-rh-parsorbitalis
60    ctx-rh-parstriangularis
61    ctx-rh-pericalcarine
62    ctx-rh-postcentral
63    ctx-rh-posteriorcingulate
64    ctx-rh-precentral
65    ctx-rh-precuneus
66    ctx-rh-rostralanteriorcingulate
67    ctx-rh-rostralmiddlefrontal
68    ctx-rh-superiorfrontal
69    ctx-rh-superiorparietal
70    ctx-rh-superiortemporal
71    ctx-rh-supramarginal
72    ctx-rh-frontalpole
73    ctx-rh-temporalpole
74    ctx-rh-transversetemporal
75    ctx-rh-insula
76    Right-Thalamus-Proper
77    Right-Caudate
78    Right-Putamen
79    Right-Pallidum
80    Right-Hippocampus
81    Right-Amygdala
82    Right-Accumbens-area
%------------------------------------------------------------------------------
Region labels for HCP parcellation can be found in Glasser et al. 2016 (doi:  10.1038/nature18933), Neuroanatomical Supplementary Results.
