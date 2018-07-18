
function figure5C()% make a histogram for assignment distances
% load any file where samples were assigned to parcellation
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC82DistThresh2.mat')

sides = {'left', 'right'};
brainParts = {'Cortex','Subcortex'};

k=1;
for i=1:6
    for side = sides
        for brainPart = brainParts
            
            if isfield(assignDistance{i}, (side{1}))
                A{k} = assignDistance{i}.(side{1}).(brainPart{1});
                k=k+1;
                
            end
        end
    end
    
end

B = vertcat(A{1}, A{2}, A{3}, A{4}, ...
    A{5}, A{6}, A{7}, A{8}, ...
    A{9} ,A{10}, A{11}, A{12}, ...
    A{13}, A{14}, A{15}, A{16});

j=1;
for l=0:20
    prop(j) = length(find(B<=l))/length(B);
    j=j+1;
end

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
% plot on large axes
figure; set(gcf,'color','w'); box off;
histogram(B, 27,'facecolor',nice_cmap(20,:),'facealpha',.5,'edgecolor',nice_cmap(10,:));
xlabel('Distance (mm)'); ylabel('Number of samples');
set(gca,'FontSize', 18)
% create smaller axes in top right, and plot on it
axes('Position',[.35 .35 .5 .5])

%box on
x=0:20;
plot(x, prop, 'Color',nice_cmap(10,:), 'LineWidth',1.5);
xlabel('Distance (mm)'); ylabel('Proportion of assigned samples');
hold on; box off;
scatter(x, prop, 'MarkerEdgeColor',nice_cmap(10,:),...
    'MarkerFaceColor',nice_cmap(20,:));
set(gca,'FontSize', 16)
cd ../../..
end

