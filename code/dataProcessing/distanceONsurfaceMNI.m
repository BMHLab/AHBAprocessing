function shortestDist = distanceONsurfaceMNI(data_parcel, ROIs, coordsAssigned)

%samplesIND = ismember(DataCoordinates{subject}(:,2),nROIs);
samples = coordsAssigned; %DataCoordinates{subject}(samplesIND,3:5);

data = ismember(data_parcel,ROIs);
ind = 1:(size(data,1)*size(data,2)*size(data,3));

dataNew = reshape(ind,size(data));
nodes = length(nonzeros(data(:)));

nodeList = zeros(nodes*3,3);

a=1;             % distance between adjacent voxels
c=sqrt(a^2+a^2); % distance between diagonal voxels according to pythagoras.
i=1;             % counting index to
% assign distance values to the selected plain
surrValues2 = [c,a,c; a,0,a; c,a,c];
% assign distance values to the more distant planes
surrValues13 = [c,c,c; c,a,c; c,c,c];

% make a 3D cube of distances
surrValues3D(1,:,:) = surrValues13;
surrValues3D(2,:,:) = surrValues2;
surrValues3D(3,:,:) = surrValues13;

% loop for each voxel in the parcellation
tic
for x=1:size(data,1)
    for y=1:size(data,2)
        for z=1:size(data,3)
            node = data(x,y,z);
            % for every nonzero node in the parcellation
            if node~=0
                % get edge weights based on distance between nodes
                weights = surrValues3D.*data(x-1:x+1,y-1:y+1,z-1:z+1);
                % select indexes for each neighbour for a list
                neighbourIDs = dataNew(x-1:x+1,y-1:y+1,z-1:z+1);
                listIDs = neighbourIDs(:);
                % assign values for each nonzero node
                nodeList(i:i+26,1) = dataNew(x,y,z);        % index of the first node
                nodeList(i:i+26,2) = listIDs;               % index of the second node
                nodeList(i:i+26,3) = weights(:);            % weight between them
                i=i+27;
            end
        end
    end
    fprintf('X loop nr %d out of %d\n', x, size(data,1))
end
toc
% remove edges with zero weight
indices = nodeList(:,3)==0;
nodeList(indices,:) = [];

% remove duplicate edges
nodeINDlist = nodeList(:,1:2); % select IDs for samples on parcellation
weightList = nodeList(:,3); % select weights for samples on parcellation

[~, loc]=ismember(nodeINDlist(:,[2 1]),nodeINDlist,'rows');
nodeINDlist = nodeINDlist(loc>=(1:size(loc,1))',:);
weightList = weightList(loc>=(1:size(loc,1))',:);
nodeListUnique(:,1:2) = nodeINDlist;
nodeListUnique(:,3) = weightList;

% make an undirected graph from a list of edges
G = graph(nodeListUnique(:,1),nodeListUnique(:,2),nodeListUnique(:,3));

% get index values for samples
indSamples = zeros(size(samples,1),1);
for j=1:size(samples,1)
    indSamples(j) = dataNew(samples(j,1), samples(j,2), samples(j,3));
end

% calculate distance between each pair of samples
tic
shortestDist = zeros(size(samples,1));
for samp1 = 1:size(samples,1)
    for samp2=samp1+1:size(samples,1)
        
        [~,shortestDist(samp1,samp2)] = shortestpath(G,indSamples(samp1),indSamples(samp2));
        
    end
end
toc
end