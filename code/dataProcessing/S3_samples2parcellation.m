function S3_samples2parcellation(options)

cd ('data/genes/processedData')

optionsSave = options;
probeSelections = options.probeSelections;
parcellations = options.parcellations;
useCUSTprobes = options.useCUSTprobes;
distanceThreshold = options.distanceThreshold;
signalThreshold = options.signalThreshold;
divideSamples = options.divideSamples;
excludeHippo = options.excludeHippocampus;
VARscale = options.VARscale;
VARperc = options.VARperc;
VARfilter = options.VARfilter;

sides = {'right', 'left'};
% separate according to brain part: cortex/subcortex
brainParts = {'Cortex', 'Subcortex'};
% loop throgh all 6 subjects
subjects = 1:6;
%------------------------------------------------------------------------------
% Separate samples according to Allen institute Structure names.
%------------------------------------------------------------------------------

if useCUSTprobes
    fprintf(1,'Using the data with CUST probes\n')
    startFileName = 'MicroarrayDataWITHcustProbesUpdatedXXX';
else
    fprintf(1,'Using the data without CUST probes\n')
    startFileName = 'MicroarrayDataProbesUpdated';
end

if signalThreshold==-1 && VARfilter
    QClabel = sprintf('noQCvar%s', VARscale);
elseif signalThreshold==-1 && ~VARfilter
    QClabel = 'noQC';
elseif signalThreshold>-1 && ~VARfilter
    QClabel = 'QC';
elseif signalThreshold>-1 && VARfilter
    QClabel = sprintf('QCvar%s', VARscale);
end



fprintf('Separating samples according to side and structure\n')
Ontology = importOntology('Ontology.csv');
regions = Ontology.structure_id_path;
cortexList = {'cortex', 'gyrus', 'gyri', 'sulcus', 'occipital pole', 'planum temporale', 'cuneus', ...
    'frontal pole', 'operculum', 'planum polare', 'temporal pole' , 'paracentral lobule' };%
% load data file
for t = 1:length(probeSelections)

    FileName = sprintf('%s%s%s.mat',startFileName, probeSelections{t}, QClabel);
    load(FileName);

    expressionSubjects = cell(6,1);
    sampleInfoSubjects = cell(6,1);

    for subject = subjects

        % select words that belong to cortex

        % sepect samples according to parts and sides

        for side = sides
            for brainPart = brainParts
                % rename variables to keep the original
                expressionS = expressionAll{subject};
                sampleInformationS = sampleInfo{subject};

                % check each structure name for side and brain pars separation
                for i=1:length(sampleInformationS.StructureNames)

                    StructureNumb = num2str(sampleInformationS.StructureID(i));
                    StructureName = sampleInformationS.StructureNames{i};

                    if strcmp(divideSamples, 'ontology') || excludeHippo
                    % create a variable that contains ontology path to be
                    % further used in filtering samples
                    indexOnt = cell(length(regions),1);
                    for j = 1:length(regions)
                        indexOnt{j} = strfind(regions{j},StructureNumb);
                    end
                    % find the last ID in the scructure path - this is the
                    % structure ID for the sample
                    fstruct = find(~cellfun(@isempty,indexOnt));
                    isStruct = zeros(length(fstruct),1);
                    for fl=1:length(fstruct)
                        idpath = regions{fstruct(fl)};
                        sl = find(idpath == '/');
                        lastP1 = sl(end)-1;
                        lastP2 = sl(length(sl)-1)+1;
                        lastID = idpath(lastP2:lastP1);
                        % check if this number matches the selected
                        % structure number
                        isStruct(fl) = strcmp(lastID,StructureNumb);
                    end
                    % get the ontology annotation to check for relevant structure IDs
                    % 4008 - cortical samples
                    % 4275 - subcortical samples
                    % 4392 - thalamus samples
                    % 4249 - hippocampal samples
                    structPath = regions{fstruct(isStruct==1)};
                    end

                    switch divideSamples
                        case 'ontology'
                            isCortex = strfind(structPath,'4008');
                            if ~isempty(isCortex)
                                corticalSample = true;
                            else
                                corticalSample = false;
                            end
                            isSubcortex = [logical(strfind(structPath,'4275')),logical(strfind(structPath,'4392'))];

                            if ~isempty(isSubcortex)
                                subcorticalSample = true;
                            else
                                subcorticalSample = false;
                            end

                        case 'listCortex'
                            index = cell(length(cortexList),1);
                            for j = 1:length(cortexList)
                                index{j} = strfind(StructureName, cortexList{j});
                            end

                            %check if cortex related name part was found
                            emptyCells = cellfun(@isempty,index);
                            if sum(emptyCells) == length(cortexList)
                                subcorticalSample = true;
                                corticalSample = false;
                            else
                                subcorticalSample = false;
                                corticalSample = true;
                            end
                    end

                    if excludeHippo
                        isHippo = strfind(structPath,'4249'); % is the sample related to hippocampus
                        hippoSample = ~isempty(isHippo);
                        if hippoSample
                            expressionS(i,:) = NaN;
                            sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                            sampleInformationS.MMCoordinates(i,:) = NaN;
                            sampleInformationS.StructureNames{i} = 'remove';
                            sampleInformationS.StructureID(i,:) = NaN;
                        end
                    end

                    % there are several samples that are not labeled 'left',
                    % 'right' - like corpus {'corpus callosum'}, so we exclude them.
                    indexLR = strfind(StructureName, side{1});

                    if strcmp(brainPart{1}, 'Cortex') % will exclude subcortical samples
                        if ~corticalSample || isempty(indexLR)  % && isempty(index2)

                            expressionS(i,:) = NaN;
                            sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                            sampleInformationS.MMCoordinates(i,:) = NaN;
                            sampleInformationS.StructureNames{i} = 'remove';
                            sampleInformationS.StructureID(i) = NaN;
                        end

                    elseif strcmp(brainPart{1}, 'Subcortex') % will exclude cortical samples

                        if  ~subcorticalSample || isempty(indexLR)
                            expressionS(i,:) = NaN;
                            sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                            sampleInformationS.MMCoordinates(i,:) = NaN;
                            sampleInformationS.StructureNames{i} = 'remove';
                            sampleInformationS.StructureID(i) = NaN;

                        end
                    end

                end

                % exclude nonexisting data
                expressionS(any(isnan(expressionS),2),:) = [];
                expression.(side{1}).(brainPart{1}) = expressionS;

                sampleInformationS.MRIvoxCoordinates(any(isnan(sampleInformationS.MRIvoxCoordinates),2),:) = [];
                sampleInformation.(side{1}).(brainPart{1}).MRIvoxCoordinates = sampleInformationS.MRIvoxCoordinates;

                sampleInformationS.MMCoordinates(any(isnan(sampleInformationS.MMCoordinates),2),:) = [];
                sampleInformation.(side{1}).(brainPart{1}).MMCoordinates = sampleInformationS.MMCoordinates ;

                sampleInformationS.StructureID(any(isnan(sampleInformationS.StructureID),2),:) = [];
                sampleInformation.(side{1}).(brainPart{1}).StructureID = sampleInformationS.StructureID;

                sampleInformationS.StructureNames(strcmp('remove',sampleInformationS.StructureNames)) = [];
                sampleInformation.(side{1}).(brainPart{1}).StructureNames = sampleInformationS.StructureNames;

            end
        end


        expressionSubjects{subject} = expression;
        sampleInfoSubjects{subject} = sampleInformation;
    end
end
cd ../../..
%------------------------------------------------------------------------------
% assign samples to parcellation
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Create variables to save data in
%------------------------------------------------------------------------------
DataExpression = cell(length(subjects),1);
DataCoordinatesMRI = cell(length(subjects),1);
DataCoordinatesMNI = cell(length(subjects),1);
%------------------------------------------------------------------------------
% Select variables according to side/brain part selections
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Do assignment for all subjects
%------------------------------------------------------------------------------
for parcellation = parcellations
    for t=1:length(probeSelections)
        assignDistance = cell(length(subjects),1);
        for subject = subjects
            cd ('data/genes/parcellations')
            subjectDir = sprintf('S0%d_H0351', subject);
            cd (subjectDir)

            if strcmp(parcellation, 'aparcaseg') || strcmp(parcellation, 'cust100') || strcmp(parcellation, 'cust250')
                brainParts = {'Cortex','Subcortex'};
            else
                brainParts = {'Cortex'};
            end
            fprintf('Subject %u parcellation %s assignment distance threshold %u\n; ', subject, parcellation{1}, distanceThreshold )

            %------------------------------------------------------------------------------
            % Load parcellations
            %------------------------------------------------------------------------------
             if strcmp(parcellation, 'aparcaseg')
                [~, data_parcel]=read('defaultparc_NativeAnat.nii');
                NumNodes = 82;
                LeftCortex = 1:34;
                LeftSubcortex = 35:41;
                RightCortex = 42:75;
                RightSubcortex = 76:82;
            elseif strcmp(parcellation, 'cust100')
                [~, data_parcel]=read('random200_acpc_uncorr_asegparc_NativeAnat.nii');
                NumNodes = 220;
                LeftCortex = 1:100;
                LeftSubcortex = 101:110;
                RightCortex = 111:210;
                RightSubcortex = 211:220;
            elseif strcmp(parcellation, 'cust250')
                [~, data_parcel]=read('random500_acpc_uncorr_asegparc_NativeAnat.nii');
                NumNodes = 530;
                LeftCortex = 1:250;
                LeftSubcortex = 251:265;
                RightCortex = 266:515;
                RightSubcortex = 516:530;
            elseif strcmp(parcellation, 'HCP')
                [~, data_parcel]=read('HCPMMP1_acpc_uncorr.nii');
                NumNodes = 360;
                LeftCortex = 1:180;
                RightCortex = 181:360;
            elseif strcmp(parcellation, 'Schaefer100')
                [~, data_parcel]=read('Schaefer100_7net_config_uncorr.nii');
                NumNodes = 100;
                LeftCortex = 1:NumNodes/2;
                RightCortex = 1:NumNodes/2+1:NumNodes;
            elseif strcmp(parcellation, 'Schaefer300')
                [~, data_parcel]=read('Schaefer300_7net_config_uncorr.nii');
                NumNodes = 300;
                LeftCortex = 1:NumNodes/2;
                RightCortex = 1:NumNodes/2+1:NumNodes;
            elseif strcmp(parcellation, 'Schaefer500')
                [~, data_parcel]=read('Schaefer500_7net_config_uncorr.nii');
                NumNodes = 500;
                LeftCortex = 1:NumNodes/2;
                RightCortex = 1:NumNodes/2+1:NumNodes;
            elseif strcmp(parcellation, 'Schaefer1000')
                [~, data_parcel]=read('Schaefer1000_7net_config_uncorr.nii');
                NumNodes = 1000;
                LeftCortex = 1:NumNodes/2;
                RightCortex = 1:NumNodes/2+1:NumNodes;
            elseif strcmp(parcellation, 'HCPmni')
                [~, data_parcel]=read('MMPinMNI.nii');
                NumNodes = 360;
                LeftCortex = 1:180;
                RightCortex = 181:360;
            end
            cd ../../

            %------------------------------------------------------------------------------
            % Load microarray data
            %------------------------------------------------------------------------------
            cd ('processedData');

            if useCUSTprobes
                fprintf(1,'Using the data with CUST probes %s for %d\n', probeSelections{t}, subject)
            else
                fprintf(1,'Using the data without CUST probes %s for %d\n', probeSelections{t}, subject)
            end

            for side = sides
                for brainPart = brainParts

                    coords2assign = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MRIvoxCoordinates;
                    coords2assignMNI = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MMCoordinates;

                    %------------------------------------------------------------------------------
                    % Find coordinates for nonzero elements in parcellation according to side and brain part (will be used to assign microarray samples to)
                    %------------------------------------------------------------------------------
                    % get coordinate values
                    if strcmp(side, 'left') && strcmp(brainPart, 'Cortex')
                        Text = sprintf('Subject %d LEFT cortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, LeftCortex)));
                    elseif strcmp(side, 'left') && strcmp(brainPart, 'Subcortex')
                        Text = sprintf('Subject %d LEFT subcortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, LeftSubcortex)));
                    elseif strcmp(side, 'right') && strcmp(brainPart, 'Cortex')
                        Text = sprintf('Subject %d RIGHT cortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, RightCortex)));
                    elseif strcmp(side, 'right') && strcmp(brainPart, 'Subcortex')
                        Text = sprintf('Subject %d RIGHT subcortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, RightSubcortex)));
                    end

                    coordsNonzeroParcel = cat(2,Coordx,Coordy,Coordz);

                    %------------------------------------------------------------------------------
                    % For each microarray coordinate find a closest coordinate in parcellation
                    %------------------------------------------------------------------------------
                    if strcmp(parcellation, 'HCPmni')
                        for samp=1:length(coords2assign)
                            coords2assign(samp,:) = mni2orFROMxyz(coords2assignMNI(samp,:),1,'mni');
                        end
                    end
                    coordsAssigned = zeros(length(coords2assign),3);
                    % find closest point
                    k = dsearchn(coordsNonzeroParcel,coords2assign);

                    for i = 1:size(k,1)
                        coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
                    end

                    %------------------------------------------------------------------------------
                    % Calculate the distance between original and reassigned coordinate
                    %------------------------------------------------------------------------------

                    for j=1:size(coordsAssigned,1)
                        assignDistance{subject}.(side{1}).(brainPart{1})(j,1) = pdist2(coordsAssigned(j,:), coords2assign(j,:));
                        if assignDistance{subject}.(side{1}).(brainPart{1})(j,1)>distanceThreshold
                            coordsAssigned(j,:)=NaN;
                        end
                    end

                    %------------------------------------------------------------------------------
                    % Extract intensity values for assigned coordinates
                    %------------------------------------------------------------------------------

                    idxAssigned = find(~any(isnan(coordsAssigned), 2));
                    intensity_all=zeros(length(idxAssigned), 1);

                    for i=1:length(idxAssigned)

                        sampIDX = idxAssigned(i);
                        intensity = data_parcel(coordsAssigned(sampIDX,1), coordsAssigned(sampIDX,2), coordsAssigned(sampIDX,3));
                        intensity_all(i) = intensity;
                    end

                    MRIcoordinates = coordsAssigned;
                    MNIcoordinates = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MMCoordinates;
                    expr = expressionSubjects{subject}.(side{1}).(brainPart{1});

                    MNIcoordinates = MNIcoordinates(idxAssigned, :);
                    MRIcoordinates = MRIcoordinates(idxAssigned, :);
                    expr = expr(idxAssigned,:);
                    %------------------------------------------------------------------------------
                    % sort nonzero coordinates according to intensity value
                    %------------------------------------------------------------------------------
                    informationMRIcoord = [intensity_all MRIcoordinates];
                    informationMNIcoord = [intensity_all MNIcoordinates];
                    informationExpression = [intensity_all expr];

                    informationMRIcoord(any(isnan(informationMRIcoord),2),:) = [];
                    informationMNIcoord(any(isnan(informationMNIcoord),2),:) = [];
                    informationExpression(any(isnan(informationExpression),2),:) = [];

                    % sort samples according to ROIs
                    [~,indSORT] = sortrows(informationMRIcoord,1);
                    informationMRIcoord = informationMRIcoord(indSORT,:);
                    informationMNIcoord = informationMNIcoord(indSORT,:);
                    informationExpression = informationExpression(indSORT,:);

                    fprintf(1,'Sorts expression data according to ROIs\n')

                    data.(side{1}).(brainPart{1}).expression = informationExpression;
                    data.(side{1}).(brainPart{1}).informationMRI = informationMRIcoord;
                    data.(side{1}).(brainPart{1}).informationMNI = informationMNIcoord;

                end

            end


            %------------------------------------------------------------------------------
            % Save output
            %------------------------------------------------------------------------------
            if strcmp(parcellation, 'aparcaseg') || strcmp(parcellation, 'cust100') || strcmp(parcellation, 'cust250')
                nSamples = size(data.left.Cortex.informationMRI,1)+size(data.left.Subcortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1)+size(data.right.Subcortex.informationMRI,1);
                Expression = cat(1,data.left.Cortex.expression,data.left.Subcortex.expression, data.right.Cortex.expression, data.right.Subcortex.expression);
                CoordinatesMRI = cat(1,data.left.Cortex.informationMRI,data.left.Subcortex.informationMRI, data.right.Cortex.informationMRI, data.right.Subcortex.informationMRI);
                CoordinatesMNI = cat(1,data.left.Cortex.informationMNI,data.left.Subcortex.informationMNI, data.right.Cortex.informationMNI, data.right.Subcortex.informationMNI);
            else
                nSamples = size(data.left.Cortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1);
                Expression = cat(1,data.left.Cortex.expression, data.right.Cortex.expression);
                CoordinatesMRI = cat(1,data.left.Cortex.informationMRI, data.right.Cortex.informationMRI);
                CoordinatesMNI = cat(1,data.left.Cortex.informationMNI, data.right.Cortex.informationMNI);
            end
            SUBJECT = zeros(nSamples,1);
            SUBJECT(:,1) = subject;

            DataExpression{subject} = [SUBJECT, Expression];
            DataCoordinatesMRI{subject} = [SUBJECT, CoordinatesMRI];
            DataCoordinatesMNI{subject} = [SUBJECT, CoordinatesMNI];

            cd ../../..
        end
        options = optionsSave;

        %% save data for all subjects
        cd ('data/genes/processedData')

        save(sprintf('%s%s%s%dDistThresh%d.mat', startFileName, probeSelections{t}, QClabel, NumNodes, distanceThreshold), 'DataExpression', 'DataCoordinatesMRI', 'DataCoordinatesMNI', 'probeInformation', 'assignDistance', 'options');

        cd ../../..

    end
end

end
