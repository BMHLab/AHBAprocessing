
%% This script:
%   1. Chooses one vertex per sample on a surface
%   2. calculates distances between samples on the surface
%%
%------------------------------------------------------------------------------
% Choose single vertex per sample based on the max value
%------------------------------------------------------------------------------

cd ('data/genes/forFreesurferXXX')
% load flsaverage brain surface
[verticesFSaverage,facesFSaverage] = read_surf('lhfsaverage.pial');
listChosenVert = cell(6,1); 

% for each subject and each sample choose a single vertex per sample
for sub=1:6
    fprintf('Processing subject %d\n', sub); 
    cd (sprintf('S0%d', sub))
    samples = load(sprintf('S%dsampleListXXX.txt', sub));
    subjChosenVert = zeros(length(samples),1);
    k=1;
    for i=1:length(samples)
    sample = samples(i);

    dataOrig = MRIread(sprintf('S%dsample%dONfsaveragePialXXX.mgz',sub, sample));
    [~,chosenVert] = max(dataOrig.vol);
    subjChosenVert(k) = chosenVert;
    k=k+1;
    end
    listChosenVert{sub} = subjChosenVert; 
    cd ..

end

% make a combined list of vertices
vertices = vertcat(listChosenVert{1}, listChosenVert{2},listChosenVert{3}, listChosenVert{4}, listChosenVert{5}, listChosenVert{6}); 
coordinates = verticesFSaverage(vertices,:); 

%------------------------------------------------------------------------------
% Calculate distances on the surface for each sample 
%------------------------------------------------------------------------------
% takes about an hour
distances = zeros(size(verticesFSaverage,1),length(vertices)); 
tic
for samp=1:length(vertices)
    
    [D,S,Q] = perform_fast_marching_mesh(verticesFSaverage.', facesFSaverage.' + 1, vertices(samp));
    distances(:,samp) = D; 
    fprintf('%d sample has been analised\n', samp)
    
end
toc
% select distances between relevant samples
distSamples = distances(vertices,:); 
cd ..
cd 'processedData'
save('DistancesONsurfaceXXX.mat', 'distSamples')
% plot distances for a random sample - produces a nice picture
figure;patch('Vertices',verticesFSaverage,'faces',facesFSaverage+1,'CData',distances(:,500),'FaceColor','interp','EdgeColor','none');axis off;camlight;axis image;view([0 45])

