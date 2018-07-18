function S3_samples2parcellation(options)

cd ('data/genes/processedData')


probeSelections = options.probeSelections;
parcellations = options.parcellations;
useCUSTprobes = options.useCUSTprobes;
distanceThreshold = options.distanceThreshold;
signalThreshold = options.signalThreshold;

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

if signalThreshold==-1
    QClabel = 'noQC'; 
else
    QClabel = 'QC'; 
end


fprintf('Separating samples according to side and structure\n')
% load data file
for t = 1:length(probeSelections)
    
    FileName = sprintf('%s%s%s.mat',startFileName, probeSelections{t}, QClabel);
    load(FileName);
    
    expressionSubjects = cell(6,1);
    sampleInfoSubjects = cell(6,1);
    
    for subject = subjects
        
        % select words that belong to cortex
        cortex = {'cortex', 'gyrus', 'gyri', 'sulcus', 'occipital pole', 'planum temporale', 'cuneus', ...
            'frontal pole', 'operculum', 'planum polare', 'temporal pole' , 'paracentral lobule' };%
        
        % sepect samples according to parts and sides
        index = cell(length(cortex),1);
        for side = sides
            for brainPart = brainParts
                % rename variables to keep the original
                expressionS = expressionAll{subject};
                sampleInformationS = sampleInfo{subject};
                
                % chack each structure name for side and brain pars separation
                for i=1:length(sampleInformationS.StructureNames)
                    Structure = sampleInformationS.StructureNames{i};
                    
                    %check if structure name has a part related to cortex
                    for j = 1:length(cortex)
                        index{j} = strfind(Structure, cortex{j});
                    end
                    
                    %check if cortex related name part was found
                    emptyCells = cellfun(@isempty,index);
                    
                    if strcmp(brainPart{1}, 'Cortex') % will exclude subcortical samples
                        
                        if sum(emptyCells) == length(cortex)% && isempty(index2)
                            expressionS(i,:) = NaN;
                            sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                            sampleInformationS.MMCoordinates(i,:) = NaN;
                            sampleInformationS.StructureNames{i} = 'remove';
                        end
                        
                    elseif strcmp(brainPart{1}, 'Subcortex') % will exclude cortical samples
                        
                        if sum(emptyCells) < length(cortex)% && isempty(index2)
                            expressionS(i,:) = NaN;
                            sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                            sampleInformationS.MMCoordinates(i,:) = NaN;
                            sampleInformationS.StructureNames{i} = 'remove';
                            
                        end
                    end
                    % there are several samples that are not labeled 'left',
                    % 'right' - like corpus {'corpus callosum'}, so we exclude them.
                    index2 = strfind(Structure, side{1});
                    if isempty(index2)
                        expressionS(i,:) = NaN;
                        sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                        sampleInformationS.MMCoordinates(i,:) = NaN;
                        sampleInformationS.StructureNames{i} = 'remove';
                    end
                end
                
                % exclude nonexisting data
                expressionS(any(isnan(expressionS),2),:) = [];
                expression.(side{1}).(brainPart{1}) = expressionS;
                
                sampleInformationS.MRIvoxCoordinates(any(isnan(sampleInformationS.MRIvoxCoordinates),2),:) = [];
                sampleInformation.(side{1}).(brainPart{1}).MRIvoxCoordinates = sampleInformationS.MRIvoxCoordinates;
                
                sampleInformationS.MMCoordinates(any(isnan(sampleInformationS.MMCoordinates),2),:) = [];
                sampleInformation.(side{1}).(brainPart{1}).MMCoordinates = sampleInformationS.MMCoordinates ;
                
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
        coordsAssignedALL = cell(length(subjects),1);
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
                    
                    if strcmp(parcellation, 'HCP')
                        brainParts = {'Cortex'};
                    else
                        brainParts = {'Cortex','Subcortex'};
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
                        [~, data_parcel]=read('customparc100_NativeAnatFixed.nii');
                        NumNodes = 220;
                        LeftCortex = 1:100;
                        LeftSubcortex = 101:110;
                        RightCortex = 111:210;
                        RightSubcortex = 211:220;
                    elseif strcmp(parcellation, 'cust250')
                        [~, data_parcel]=read('customparc250_NativeAnatFixed.nii');
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
                            
                            if distanceThreshold < 35
                                
                                load(sprintf('CoordsAssignedAll%d.mat', NumNodes));
                                
                            end
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
                            
                            coordsAssigned =  zeros(size(coords2assign));
                            % find closest point
                            k = dsearchn(coordsNonzeroParcel,coords2assign);
                            
                            for i = 1:size(k,1)
                                coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
                            end
                            
                            %------------------------------------------------------------------------------
                            % Salculate the distance between original and reassigned coordinate
                            %------------------------------------------------------------------------------
                            %(size(coordsAssigned,1),1);
                            coordsNONassigned = zeros(size(coordsAssigned,1),3);
                            coordsToAssignAll = coords2assign;
                          
                            for j=1:size(coordsAssigned,1)
                                
                                assignDistance{subject}.(side{1}).(brainPart{1})(j,1) = pdist2(coordsAssigned(j,:), coords2assign(j,:));
                                
                                if assignDistance{subject}.(side{1}).(brainPart{1})(j,1)>distanceThreshold
                                    coordsNONassigned(j,:,:,:)=coords2assign(j,:,:,:);
                                    coordsAssigned(j,:)=0;
                                    coords2assign(j,:)=0;
                                    
                                end
                                
                            end
                            %remove zero elements
                            coords2assign( ~any(coords2assign,2), : ) = NaN;
                            coordsAssigned( ~any(coordsAssigned,2), : ) = NaN;
                            
                            coordsNONassigned( ~any(coordsNONassigned,2), : ) = [];
                            
                            %------------------------------------------------------------------------------
                            % Extract intensity values for assigned coordinates
                            %------------------------------------------------------------------------------
                            
                            if distanceThreshold >35
                                
                                coordsAssignedALL{subject}.(side{1}).(brainPart{1}) = coordsAssigned;
                                
                            elseif distanceThreshold <35
                                
                                Int=nonzeros(data_parcel);
                                int = unique(Int);
                                
                                intensity_all=zeros(size(coordsAssignedALL{subject}.(side{1}).(brainPart{1}),1), 1);
                                
                                for i=1:size(coordsAssignedALL{subject}.(side{1}).(brainPart{1}),1)
                                    
                                    intensity = data_parcel(coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,1), coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,2), coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,3));
                                    intensity_all(i) = intensity;
                                    
                                end
                                indRemove = find(any(isnan(coordsAssigned),2));
                                intensity_all(indRemove, :) = [];
                                
                                MNIcoordinates = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MMCoordinates;
                                
                                MNIcoordinates(indRemove, :) = [];
                                coordsAssigned(indRemove, :) = [];
                                
                                %------------------------------------------------------------------------------
                                % sort nonzero coordinates according to intensity value
                                %------------------------------------------------------------------------------
                                informationALL = [intensity_all coordsAssigned];
                                
                                
                                informationALLmni = [intensity_all MNIcoordinates];
                                information_all = informationALL;
                                informationmni = informationALLmni;
                                
                                informationALL(any(isnan(informationALL),2),:) = NaN;
                                information_all(any(isnan(informationALL),2),:) = [];
                                
                                informationALLmni(any(isnan(informationALL),2),:) = NaN;
                                informationmni(any(isnan(informationALL),2),:) = [];
                                
                                information_sorted = sortrows(information_all,1);
                                information_sortedmni = sortrows(informationmni,1);
                                
                                
                                sampleIntensity = informationALL(:,1);
                                expr = expressionSubjects{subject}.(side{1}).(brainPart{1});
                                expr(indRemove,:)=[];
                                IntANDExpression = [sampleIntensity expr];
                                % exclude nonassigned points
                                IntANDExpression(any(isnan(IntANDExpression),2),:) = NaN;
                                IntANDExpression(any(isnan(IntANDExpression),2),:) = [];
                                
                                fprintf(1,'Sorts expression data according to ROIs\n')
                                
                                % sort according to intensity (roi)
                                IntANDExpressionSorted = sortrows(IntANDExpression,1);
                                % assign outputs to structure
                                information.(side{1}).(brainPart{1}).raw = information_all;
                                information.(side{1}).(brainPart{1}).sorted = information_sorted;
                                coordinates.(side{1}).(brainPart{1}).assigned = coordsAssigned;
                                coordinates.(side{1}).(brainPart{1}).toassign = coords2assign;
                                coordinates.(side{1}).(brainPart{1}).toassignALL = coordsToAssignAll;
                                coordinates.(side{1}).(brainPart{1}).nonassigned = coordsNONassigned;
                                
                                data.(side{1}).(brainPart{1}).expression = IntANDExpressionSorted;
                                data.(side{1}).(brainPart{1}).informationMRI = information_sorted;
                                data.(side{1}).(brainPart{1}).informationMNI = information_sortedmni;
                                
                            end
                            
                            %% plot
                            %     figure; scatter3(coordsAssigned(:,1), coordsAssigned(:,2), coordsAssigned(:,3), 'bo');...
                            %         hold on; scatter3(coordsToAssign(:,1), coordsToAssign(:,2), coordsToAssign(:,3), 'g*');...
                            %         hold on; scatter3(coordsNonassigned(:,1), coordsNonassigned(:,2), coordsNonassigned(:,3), 'r.');
                            %         hold on;
                            %         % plot arrows where the point moved
                            %         for j = 1:length(coordsToAssign)
                            %             dp = coordsToAssign(j,:) - coordsAssigned(j,:);
                            %             q = quiver3(coordsToAssign(j,1),coordsToAssign(j,2),coordsToAssign(j,3),-dp(1), -dp(2), -dp(3), 'LineWidth', 1.5);
                            %             hold on;
                            %         end
                            
                            %savefig(sprintf('%d_DistThresh%d_coordMapping_MRI.fig', NumNodes, DistanceThreshold ));
                            
                            %close all;
                        end
                        
                    end
                    
                    
                    %------------------------------------------------------------------------------
                    % Save output
                    %------------------------------------------------------------------------------
                    if distanceThreshold <35
                        if strcmp(parcellation, 'HCP')
                            nSamples = size(data.left.Cortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1);
                            Expression = cat(1,data.left.Cortex.expression, data.right.Cortex.expression);
                            CoordinatesMRI = cat(1,data.left.Cortex.informationMRI, data.right.Cortex.informationMRI);
                            CoordinatesMNI = cat(1,data.left.Cortex.informationMNI, data.right.Cortex.informationMNI);
                        else
                            nSamples = size(data.left.Cortex.informationMRI,1)+size(data.left.Subcortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1)+size(data.right.Subcortex.informationMRI,1);
                            Expression = cat(1,data.left.Cortex.expression,data.left.Subcortex.expression, data.right.Cortex.expression, data.right.Subcortex.expression);
                            CoordinatesMRI = cat(1,data.left.Cortex.informationMRI,data.left.Subcortex.informationMRI, data.right.Cortex.informationMRI, data.right.Subcortex.informationMRI);
                            CoordinatesMNI = cat(1,data.left.Cortex.informationMNI,data.left.Subcortex.informationMNI, data.right.Cortex.informationMNI, data.right.Subcortex.informationMNI);
                        end
                        SUBJECT = zeros(nSamples,1);
                        SUBJECT(:,1) = subject;
                        
                        DataExpression{subject} = [SUBJECT, Expression];
                        DataCoordinatesMRI{subject} = [SUBJECT, CoordinatesMRI];
                        DataCoordinatesMNI{subject} = [SUBJECT, CoordinatesMNI];
 
                    elseif distanceThreshold > 35
                        
                        save(sprintf('CoordsAssignedAll%d.mat', NumNodes), 'coordsAssignedALL');
                        
                    end
                    cd ../../..
                end
                
                
                %% save data for all subjects
                cd ('data/genes/processedData')
                if distanceThreshold < 35
                    save(sprintf('%s%s%s%dDistThresh%d.mat', startFileName, probeSelections{t}, QClabel, NumNodes, distanceThreshold), 'DataExpression', 'DataCoordinatesMRI', 'DataCoordinatesMNI', 'probeInformation', 'assignDistance', 'options');
                end
                cd ../../..
                
            end
        end
        
end
