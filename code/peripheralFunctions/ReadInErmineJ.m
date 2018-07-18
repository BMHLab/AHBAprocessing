function [summaryTable, GOName,GOID,pval,corr_pval,numGenes,geneMembers] = ReadInErmineJ(fileNameIn,fileNameMatch)
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-12-10
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
if nargin < 2
    matchFileName = 0;
else
    matchFileName = 1;
end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Open the file
% ------------------------------------------------------------------------------
% fileName = 'connected_coor_tstat.project');
fid = fopen(fileNameIn,'r');

if matchFileName
    % Start scanning for the matching scorefile:
    stopHere = 'scoreFile';
    stopNow = 0;
    while (stopNow==0)
        C = textscan(fid,'%s',1,'Delimiter','');
        % disp(C{1}{1}(1:min(length(stopHere),length(C{1}{1}))))
        % pause(0.1)
        if strcmp(C{1}{1}(1:min(length(stopHere),length(C{1}{1}))),stopHere);
            % C = textscan(fid,'%s',3);
            % fileNameThisAnal = C{1};
            % Does it match?
            parts = regexp(C{1},'/','split');
            theFileName = parts{1}{end};
            if strcmp(theFileName,fileNameMatch)
                fprintf(1,'Match to %s\n',theFileName);
                % This is it!!
                stopNow = 1;
            else
                fprintf(1,'No match to %s\n',theFileName);
            end
        end
        if feof(fid)
            error('%s not found in %s',fileNameMatch,fileNameIn);
        end
    end
end

% So we found our match, now start scanning for the following
stopHere = '#!';
stopNow = 0;
while (stopNow==0)
    C = textscan(fid,'%s',1);
    if strcmp(C{1},stopHere);
        stopNow = 1;
    end
    if feof(fid)
        error('end of file');
    end
end

% ------------------------------------------------------------------------------
% Now read headings:
% ------------------------------------------------------------------------------
% #!    Name    ID    NumProbes    NumGenes    RawScore    Pval    CorrectedPvalue    MFPvalue    CorrectedMFPvalue    Multifunctionality    Same as    GeneMembers
headings = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1);

% Now read data to the end of this block:
data1 = textscan(fid,'%s %s %s %u %u %f %f %f %f %f %f %s %s','Delimiter','\t'); %,'MultipleDelimsAsOne',1);
% !    regulation of neuronal synaptic plasticity    GO:0048168    39    37    0.73788422    1.192E-07    5.207E-04    3.342E-07    1.46E-03    0.866
fclose(fid)

fprintf(1,'Found %u rows of results in %s\n',length(data1{1}),fileNameIn);

GOName = data1{2};
GOID = data1{3};
numGenes = data1{5};
% GOIDnum = cellfun(@(x)x(4:end),GOIDstr,'UniformOutput',0);

pval = data1{7};
corr_pval = data1{8};
geneMembers = data1{13};

summaryTable = table(GOID,GOName,numGenes,pval,corr_pval,geneMembers);

% Nonsense, 'Name'    'ID'    'NumProbes'    'NumGenes'    'RawScore'    'Pval'    'CorrectedPvalue'    'MFPvalue'  'CorrectedMFPvalue'    'Multifunctionality'    'Same as GeneMembers'

end
