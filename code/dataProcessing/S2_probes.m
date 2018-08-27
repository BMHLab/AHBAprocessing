%% ------------------------------------------------------------------------------
% Save relevant variables to a MicroarrayData.mat file
%%------------------------------------------------------------------------------
function [expressionAll, probeInformation] = S2_probes(options)
% percentage of samples that a selected probe has expression levels that are higher than background
rng shuffle % for selecting different seed for random probe selection
useCUSTprobes = options.useCUSTprobes;
probeSelections = options.probeSelections{1};
signalThreshold = options.signalThreshold;
saveOutput = options.saveOutput;

if strcmp(probeSelections, 'RNAseq')
    RNAseqThreshold = options.RNAseqThreshold;
    RNAsign = options.RNAsignThreshold; 
end

if signalThreshold==-1
    QClabel = 'noQC'; 
else
    QClabel = 'QC'; 
end

%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');

if useCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    
    startFileName = 'MicroarrayDataWITHcustProbesUpdatedXXX';
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    startFileName = 'MicroarrayDataProbesUpdatedXXX';
end
fprintf(1,sprintf('Probe selection based on %s is chosen\n', probeSelections))
load(sprintf('%s.mat', startFileName));

cd ..
cd ('rawData');

%------------------------------------------------------------------------------
% Calculate probe selection criteria for each subject separately (non
% normalised data)
% Choose probe with max probe selection criteria on average
%------------------------------------------------------------------------------

%% find repeating entrezIDs and calculate variances for them, then take a probe of max variance.
%------------------------------------------------------------------------------
% Find best representative in a set of duplicates using maxVar and remove all others:
%------------------------------------------------------------------------------
expressionSelected = cell(6,1);
noiseSUBJ = cell(6,1);

ProbeID = DataTableProbe.ProbeID{1,1};
% % ------------------------------------------------------------------------------
% First, find probes that have very noisy data and remove them from consideration
% Threshold for removing those proges is defined as the percentage of
% samples a probe has expression higher that background
% % ------------------------------------------------------------------------------

signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=signalThreshold);

% remove selected probes from data and perform other calculations only on
% non-noisy probes
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
ProbeID = ProbeID(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);
GeneSymbol = DataTableProbe.GeneSymbol{1,1}(indKeepProbes);

probeInformationALL.ProbeName = ProbeName;
probeInformationALL.ProbeID = ProbeID;
probeInformationALL.EntrezID = EntrezID;
probeInformationALL.GeneSymbol = GeneSymbol;


% if choosing probes based on RNAseq, then use data only from 1 subject
if strcmp(probeSelections, 'RNAseq')
    
    [correlations, avgCorr, avgPval, indProbe, genes, overlapStructures] = selectProbeRNAseq(DataTable, EntrezID, indKeepProbes, RNAseqThreshold, RNAsign);
    nSub = 1;
elseif strcmp(probeSelections, 'DS')
    [indProbe, avCorr] = selectProbeDS(EntrezID, DataTable, indKeepProbes);
    nSub = 1;
else
    nSub = 6;
end


Uniq = unique(EntrezID);
ProbeList = zeros(length(Uniq),2);
indMsubj = zeros(length(Uniq),nSub);

for subj = 1:nSub
    expression = (DataTable.Expression{subj}(indKeepProbes,:))';
    noise = DataTable.Noise{subj}(indKeepProbes,:)';
    
    for k=1:length(Uniq)
        %fprintf(1,'Processing entrez ID %u\n',Uniq(k))
        % find indexes for repeating entrexIDs
        indRepEntrezIDs = find(EntrezID==(Uniq(k)));
        expRepEntrezIDs = expression(:,indRepEntrezIDs);
        if length(indRepEntrezIDs) >=2
            
            %fprintf(1,'%d duplicates found\n', length(indRepEntrezIDs));
            
            % take expression values for a selected entrezID
            % calculate variances for expression data for a selected entrezID
            %switch probeSelection
            if strcmp(probeSelections, 'Variance')
                %case 'Variance'
                %fprintf(1,'Performing probe selection using Max variance\n');
                measure = var(expRepEntrezIDs,0,1);
                % determine max var value
                [~, indMaxV] = max(measure);
                
                
            elseif strcmp(probeSelections, 'PC')
                %fprintf(1,'Performing probe selection using max PC\n');
                % substract the mean before doing pca
                expRepEntrezIDsNOmean = expRepEntrezIDs-mean(expRepEntrezIDs);
                measure = pca(expRepEntrezIDsNOmean,'Centered',false);
                % determine max PC loading
                [~, indMaxV] = max(measure(:,1));
            elseif strcmp(probeSelections, 'maxIntensity')
                
                measure = mean(expRepEntrezIDs);
                % determine probe with max signal
                [~, indMaxV] = max(measure);
                
            elseif strcmp(probeSelections, 'maxCorrelation_intensity')
                % accordin to the description in WGCNA, if there are 2
                % probes, one with max variance is chosen
                if size(expRepEntrezIDs,2)==2
                    measure = mean(expRepEntrezIDs);
                    % determine probe with max signal
                    [~, indMaxV] = max(measure);
                    % if there are more then 2 probes, correelations between
                    % each pair of probes is chosen and probe with max
                    % correlaition to others is selected
                elseif size(expRepEntrezIDs,2)>2
                    % calculate correlations between probes (all pairs)
                    rPR = zeros(size(expRepEntrezIDs,2),size(expRepEntrezIDs,2));
                    for pr1=1:size(expRepEntrezIDs,2)
                        for pr2=pr1+1:size(expRepEntrezIDs,2)
                            rPR(pr1, pr2) = corr(expRepEntrezIDs(:,pr1), expRepEntrezIDs(:,pr2), 'Type','Pearson');
                        end
                    end
                    rPR = rPR+rPR';
                    rPR(logical(eye(size(rPR)))) = 1;
                    A = (0.5+0.5*rPR);
                    % sum per column (or row) to get "total" correlaion
                    measure = sum(A,1);
                    [~, indMaxV] = max(measure);
                end
                
            elseif strcmp(probeSelections, 'maxCorrelation_variance')
                % accordin to the description in WGCNA, if there are 2
                % probes, one with max variance is chosen
                if size(expRepEntrezIDs,2)==2
                    measure = var(expRepEntrezIDs,0,1);
                    % determine probe with max signal
                    [~, indMaxV] = max(measure);
                    % if there are more then 2 probes, correelations between
                    % each pair of probes is chosen and probe with max
                    % correlaition to others is selected
                elseif size(expRepEntrezIDs,2)>2
                    % calculate correlations between probes (all pairs)
                    rPR = zeros(size(expRepEntrezIDs,2),size(expRepEntrezIDs,2));
                    for pr1=1:size(expRepEntrezIDs,2)
                        for pr2=pr1+1:size(expRepEntrezIDs,2)
                            rPR(pr1, pr2) = corr(expRepEntrezIDs(:,pr1), expRepEntrezIDs(:,pr2), 'Type','Pearson');
                        end
                    end
                    rPR = rPR+rPR';
                    rPR(logical(eye(size(rPR)))) = 1;
                    A = (0.5+0.5*rPR);
                    % sum per column (or row) to get "total" correlaion
                    measure = sum(A,1);
                    [~, indMaxV] = max(measure);
                end
                
            elseif strcmp(probeSelections, 'CV')
                %fprintf(1,'Performing probe selection using max PC\n');
                measure = std(expRepEntrezIDs)./mean(expRepEntrezIDs);
                % determine max PC loading
                [~, indMaxV] = max(measure);
                
            elseif strcmp(probeSelections, 'LessNoise')
                %fprintf(1,'Performing probe selection using less noise criteria\n');
                noiseRepEntrezIDs = noise(:,indRepEntrezIDs);
                % find probe with most signal in it compared to noise
                measure = sum(noiseRepEntrezIDs,1);
                % determine probe with max signal
                [~, indMaxV] = max(measure);
                
            elseif strcmp(probeSelections, 'Mean')
                expressionSelected{subj}(:,k) = mean(expRepEntrezIDs,2);
                [indMaxV] = randsample(1:size(expRepEntrezIDs,2),1);
                
                % choose one of probes as a place holder (at random)
                % for probeInformation (expression values are mean of
                % al probes)
                
            elseif strcmp(probeSelections, 'Random1')
                %fprintf(1,'Performing probe selection using random selection\n');
                
                % determine max var value
                [indMaxV] = randsample(1:size(expRepEntrezIDs,2),1);
            elseif strcmp(probeSelections, 'Random2')
                %fprintf(1,'Performing probe selection using random selection\n');
                
                % determine max var value
                [indMaxV] = randsample(1:size(expRepEntrezIDs,2),1);
            elseif strcmp(probeSelections, 'RNAseq')
                indMaxV = indProbe(k);
            elseif strcmp(probeSelections, 'DS')
                %indMaxV = indProbe(k);
                indMsubj(k,subj) = indProbe(k);
                
            end
            
            if (strcmp(probeSelections, 'Mean') || strcmp(probeSelections, 'Variance') || strcmp(probeSelections, 'Random1') ...
                    ||  strcmp(probeSelections, 'Random2') || strcmp(probeSelections, 'maxCorrelation_variance') ...
                    || strcmp(probeSelections, 'maxCorrelation_intensity') || strcmp(probeSelections, 'PC') ...
                    || strcmp(probeSelections, 'LessNoise') || strcmp(probeSelections, 'maxIntensity') ...
                    || strcmp(probeSelections, 'RNAseq') || strcmp(probeSelections, 'CV'))
               
                %if NaN, use NaN;
                if isnan(indMaxV)
                    indMsubj(k,subj) = NaN;
                    
                else
                    indMsubj(k,subj) = indRepEntrezIDs(indMaxV);
                end
            end
            
        else
            if strcmp(probeSelections, 'Mean')
                expressionSelected{subj}(:,k) = expRepEntrezIDs;
                indMsubj(k,subj) = indRepEntrezIDs;
                
            elseif strcmp(probeSelections, 'RNAseq')
                if isnan(indProbe(k))
                    indMsubj(k,subj) = NaN;
                else
                    indMsubj(k,subj) = indRepEntrezIDs;
                end
            else
                indMsubj(k,subj) = indRepEntrezIDs;
            end
            
        end
        
    end
end

for j=1:length(Uniq)
    
    indINlist = mode(indMsubj(j,:),2);
    if isnan(indINlist)
        ProbeList(j,1) = NaN;
        ProbeList(j,2) = NaN;
    else
        ProbeList(j,1) = ProbeID(indINlist);
        ProbeList(j,2) = EntrezID(indINlist);
    end
    
end
% %% check to exclude probes with not maximum variance or max PC and repeating entrezIDs (assign NaN value to
% % them)
[reordered, reorder_ind] = sort(EntrezID);
% reorder all values based on sorted entrezIDs, because this is the order
% of selected probes (1...n) from Unique.
EntrezID = EntrezID(reorder_ind);
ProbeID = ProbeID(reorder_ind);
GeneSymbol = GeneSymbol(reorder_ind);
ProbeName = ProbeName(reorder_ind);

[a,ind2rem] = setdiff(ProbeID, ProbeList(:,1));
ProbeID(ind2rem) = NaN;
EntrezID(ind2rem) = NaN;

EntrezID(isnan(EntrezID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
ProbeName(isnan(ProbeID)) = [];


probeInformation.EntrezID = EntrezID; %[reordered, reorder_ind] = sort(probeInformation.EntrezID);
probeInformation.GeneSymbol = GeneSymbol;
probeInformation.ProbeName = ProbeName;
ProbeID2nd = ProbeID; % assign probe ID values to a variable (will be used to filter gene expression values)
ProbeID(isnan(ProbeID)) = []; % remove redundant probes
probeInformation.ProbeID = ProbeID;

cd ..
cd ('processedData')
expressionAll = cell(6,1);
sampleInfo = cell(6,1);
for subject=1:6
    if strcmp(probeSelections, 'Variance') || strcmp(probeSelections, 'PC') || strcmp(probeSelections, 'maxIntensity')...
            || strcmp(probeSelections, 'LessNoise') || strcmp(probeSelections, 'Random1') || strcmp(probeSelections, 'Random2') ...
            || strcmp(probeSelections, 'maxCorrelation_intensity') || strcmp(probeSelections, 'maxCorrelation_variance') ...
            || strcmp(probeSelections, 'RNAseq') || strcmp(probeSelections, 'DS') || strcmp(probeSelections, 'CV')
        
        % exclude NaN probes keeping 1 probe for 1 entrezID.
        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        Expression = DataTable.Expression{subject,1}(indKeepProbes,:); % filter noisy probes
        Expression = Expression(reorder_ind,:); % reorder probes according to sorted entrezIDs
        Expression(isnan(ProbeID2nd),:) = []; % exclude probes that were not selected
        expressionAll{subject} = Expression';
        
    elseif strcmp(probeSelections, 'Mean')

        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        expressionAll{subject} = expressionSelected{subject};
    end
    % combine sample information variables to a structure.
    SampleInformation.StructureNames = DataTable.StructureName{subject,1};
    SampleInformation.MMCoordinates = DataTable.MMcoordinates{subject,1};
    SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};
    sampleInfo{subject} = SampleInformation;
    
end

if saveOutput
    if strcmp(probeSelections, 'RNAseq')
        
            save(sprintf('%s%s%s.mat', startFileName, probeSelections, QClabel), 'expressionAll', 'probeInformation' , 'sampleInfo', 'avgCorr', 'avgPval', 'probeInformationALL', 'genes', 'options');

    else
            save(sprintf('%s%s%s.mat', startFileName, probeSelections, QClabel), 'expressionAll', 'probeInformation' , 'sampleInfo', 'options');
    end
end
cd ../../..
