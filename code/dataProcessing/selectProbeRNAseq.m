% function to select probe based on RNA seq values
function [correlations, avgCorr, avgPval, indProbe, genes, overlapStructures] = selectProbeRNAseq(DataTable, EntrezID, indKeepProbes, threshold, thresholdSign)

numGenes = length(unique(EntrezID));
genes = unique(EntrezID);
correlations = cell(numGenes,2);
pvalues = cell(numGenes,2);
 
for subject=1:2
    
    folderName = sprintf('rnaseq_donor0%d', subject);
    cd (folderName)
    
    RNAseqinfo = importRNAseqInfo('Genes.csv');
    RNAseqTPM = importRNAseq('RNAseqTPM.csv');
    
    % replace gene names with IDs in the expression file
    RNAseqTPM = cell2mat(RNAseqTPM(:,2:end));
    RNAseqGene = RNAseqinfo.NCBIgeneID;
    
    % import info about RNAseq tissue samples
    %structureIDrna = importRNAseqStructInf('SampleAnnot.csv');
    % import information about microarray samples
    %structureIDmic = DataTable.SampleID{subject};
    % use well ID instead of SampleID (well is higher resolution)
    structureIDmic = DataTable.WellID{subject};
    structureIDrna = importRNAseqWellID('SampleAnnot.csv');
    
    % import expression values for microarray
    microarray = DataTable.Expression{subject}(indKeepProbes,:);
    microarrayGene = EntrezID;
    uniqueGenes = unique(EntrezID);
    
    % find overlapping structures between microarray and RNAseq
    overlapStructures = intersect(structureIDmic, structureIDrna);
    % average if there are multiple samples in the same structure
    expmic = zeros(size(microarray,1), length(overlapStructures));
    exprna = zeros(size(RNAseqTPM,1), length(overlapStructures));
    for struct=1:length(overlapStructures)
        % find ind for microarray and average expression for those elements
        indmic = structureIDmic==overlapStructures(struct);
        expmic(:,struct) = mean(microarray(:,indmic),2);
        % find ind for RNAseq
        indrna = structureIDrna==overlapStructures(struct);
        exprna(:,struct) = mean(RNAseqTPM(:,indrna),2);
    end
    
    % correlate each microarray probe with RNAseq of the corresponding gene, if
    % gene is not available put NaNs.
    
    for g=1:numGenes
        indRNA = find(RNAseqGene==uniqueGenes(g));
        indMIC = find(microarrayGene==uniqueGenes(g));
        % for each probe found in microarray correlate it with RNAseq
        if isempty(indRNA)
            correlations{g,subject} = NaN;
            pvalues{g,subject} = NaN;
        else
            C = zeros(length(indMIC),1);
            Pvals = zeros(length(indMIC),1);
            for p=1:length(indMIC)
                [C(p),Pvals(p)] = corr(exprna(indRNA,:)',expmic(indMIC(p),:)', 'type', 'Spearman');
               
            end
            
            correlations{g, subject} = C;
            pvalues{g, subject} = Pvals; 
            
        end
    end
    cd ..
end

% take an average of results between two subjects:
% if multiple values are available, choose the one which is on average
% higher
% keep the probe if on average between two subjects correlation is higher
% than 0.1 (ir some other threshold)
avgCorr = cell(numGenes,1);
avgPval = cell(numGenes,1);

indProbe = zeros(numGenes,1);
for gene = 1:numGenes
    cors = [correlations{gene,1}, correlations{gene,2}];
    vals = [pvalues{gene,1}, pvalues{gene,2}]; 
    avc = mean(cors,2);
    avp = mean(vals,2); 
    
    avgCorr{gene} = avc;
    avgPval{gene} = avp; 
    [chosenVal, chosenInd] = max(avc);
    if thresholdSign
    if chosenVal > threshold && avp(chosenInd)<0.05
        indProbe(gene) = chosenInd;
    else
        indProbe(gene) = NaN;
    end
    else
        if chosenVal > threshold 
        indProbe(gene) = chosenInd;
    else
        indProbe(gene) = NaN;
        end
    end
        
    
end

end

