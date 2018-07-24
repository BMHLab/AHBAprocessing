# AHBAprocessing

This repository provides the code for reproducing a range of analyses involving processing Allen Human Brain Atlas gene expression data and running a range of analyses to evaluate the effect of different processing choices.

Please read our paper, and if you use this code, please cite our paper:
A. Arnatkeviciute, B.D. Fulcher, A. Fornito. A practical guide to linking brain-wide gene expression and neuroimaging data. 

Contact Aurina Arnatkeviciute by [email](mailto:aurina.arnatkeviciute@monash.edu).

### Dependencies: 
[Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/), Version 2017-15-01.
[Toolbox Fast Marching](https://au.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching), Version 1.2.0.0 

### Scripts
In order to reproduce data used in the analyses please download the contents of this repository.

### Data files
Data files required for this project are hosted on [this figshare repository](https://figshare.com/s/441295fe494375aa0c13). 
Please create a directory called data in the root directory and download data there.

### Data processing
After retrieving data add relevant paths using AddPaths.m script. Data can be generated by running processingPipleine.m script with selected options. 
This script will run 4 data processing functions: <br />
S1_extractData - extracts the data from excel files and saves it in matlab format; <br /> 
S2_probes - performsprobe selection and data filtering if required; <br />
S3_samples2parcellation - assigns samples to the regions of interest in the parcellation; <br />
S4_normalisation - performs data normalisation and calculates correlated gene expression and also provides data in the region x gene format. 

See README_processingPipeline.txt for more information. 

### Analysis
Figures in our paper can be regenerated using makeFigures.m script in the plotting folder. All scripts should be run from the root directory. Other relevant analysis scripts are in dataProcessing.
