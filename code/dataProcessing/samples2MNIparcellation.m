% match samples using MNI coordinates
% For each sample of gene expression convert MNI coordinates into voxell
% coordinates for matching.
% take samples only from left cortex and combine them into one matrix
function coordsAssigned = samples2MNIparcellation(DataCoordinatesMNI, data_parcel, ROIs)
samples = cell(6,1);
%distanceThreshold = 10; 

for i=1:6
    allsamples = DataCoordinatesMNI{i};
    LCsamples = ismember(allsamples(:,2),1:34);
    samples{i} = allsamples(LCsamples==1,3:5);
end
data = vertcat(samples{1}, samples{2}, samples{3}, samples{4}, samples{5}, samples{6});

%for each sample in the list get voxel based coordinates
xyznew = zeros(size(data,1),3);
for samp = 1:size(data,1)
    
    xyznew(samp,:) = mni2orFROMxyz(data(samp,:), 'mni');
    
end

% load parcellation in MNI space

[Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, ROIs)));
coordsNonzeroParcel = cat(2,Coordx, Coordy, Coordz);


%------------------------------------------------------------------------------
% For each microarray coordinate find a closest coordinate in parcellation
%------------------------------------------------------------------------------

coordsAssigned = zeros(size(xyznew,1),3);
% find closest point
k = dsearchn(coordsNonzeroParcel,xyznew);

for i = 1:size(k,1)
    coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
end

%------------------------------------------------------------------------------
% Salculate the distance between original and reassigned coordinate
%------------------------------------------------------------------------------
assignDistance = zeros(size(coordsAssigned,1),1);
coordsNONassigned = zeros(size(coordsAssigned,1),3);


for j=1:size(coordsAssigned,1)
    
    assignDistance(j,1) = pdist2(coordsAssigned(j,:), xyznew(j,:));
    
%     if assignDistance(j)>distanceThreshold
%         coordsNONassigned(j,:,:,:)=xyznew(j,:,:,:);
%         coordsAssigned(j,:)=0;
%         xyznew(j,:)=0;
%         
%     end

end
%remove zero elements
xyznew( ~any(xyznew,2), : ) = NaN;
coordsAssigned( ~any(coordsAssigned,2), : ) = NaN;
coordsNONassigned( ~any(coordsNONassigned,2), : ) = [];
end