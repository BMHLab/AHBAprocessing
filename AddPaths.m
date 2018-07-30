%-------------------------------------------------------------------------------
%% Add paths
%-------------------------------------------------------------------------------
fprintf(1,'Adding all subdirectories to the Matlab path...');
% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
if isfield(directories,'folder')
    paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
else
    paths = arrayfun(@(x)fullfile(pwd,directories(x).name),1:length(directories),'UniformOutput',false);
end
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end
fprintf(1,' Added.\n');
