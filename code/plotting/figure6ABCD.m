% Select only cortical samples
function explained = figure6ABCD(doNormalise,numROIs,whatNormalisation,uselimma)

if nargin < 3
    whatNormalisation = 'scaledRobustSigmoid'; 
    uselimma = false;
end
if nargin < 4
    uselimma = false;
end

cd ('data/genes/processedData');
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqQC82DistThresh2.mat');

cd ..
cd ('rawData');

selectRegions = 1:numROIs;
D = cell(6,1);
subjNr = cell(6,1);

for s=1:6
    data = DataExpression{s};
    select = ismember(data(:,2), selectRegions);
    
    cortexData = data(select==1,3:end);
    if doNormalise
        
        expData = BF_NormalizeMatrix(cortexData,whatNormalisation);
        
    else
        expData = cortexData;
    end
    D{s} = expData;
    subjNr{s} = data(select==1,1);
    R{s} = data(select==1,2);
    
end

if uselimma
    cortexData = importlimma('limmanormalisedExpression.txt');
    cortexData = cortexData';
    if doNormalise
        expression = BF_NormalizeMatrix(cortexData,whatNormalisation);
        whatNormalisation = 'limma + scaledRobustSigmoid';
    else
      
        expression = cortexData;
    end
    
else
    
    expression = vertcat(D{1}, D{2}, D{3}, D{4}, D{5}, D{6});
end

subjects = vertcat(subjNr{1}, subjNr{2}, subjNr{3}, subjNr{4}, subjNr{5}, subjNr{6});
Regs = vertcat(R{1}, R{2}, R{3}, R{4}, R{5}, R{6});
cortS = find(Regs<=34); 
subcortS = find(Regs>34); 
pCort = cortS(randperm(length(cortS)));
pSubcort = subcortS(randperm(length(subcortS)));

if ~doNormalise
    expression = zscore(expression);
     whatNormalisation = 'non-normalized data'; 
else
    expression = expression; 
   
end
[W,score,~,~,explained] = pca(expression, 'NumComponents',8);

if numROIs==34
xC = score(pCort,1); yC = score(pCort,2); zC = score(pCort,3);
CC = subjects(pCort);
elseif numROIs==41
xC = score(pCort,1); yC = score(pCort,2); zC = score(pCort,3);
CC = subjects(pCort);
xS = score(pSubcort,1); yS = score(pSubcort,2); zS = score(pSubcort,3);
CS = subjects(pSubcort);
end

S = ones(length(subjects),1)+60;

if numROIs==34
    A = zeros(size(CC,1),3);
    
    for s=1:length(CC)
        if CC(s)==1
            A(s,:) = [0 0 .5]; % dark blue
        elseif CC(s)==2
            A(s,:) = [.53 .83 .97]; % light blue
        elseif CC(s)==3
            A(s,:) = [.7 .7 .7]; % grey
        elseif CC(s)==4
            A(s,:) = [1 .60 .40]; % orange
        elseif CC(s)==5
            A(s,:) = [.81 .07 .15]; % red
        elseif CC(s)==6
            A(s,:) = [.99 .87 .09]; % yellow
        end
    end
    
elseif numROIs==41
    C = vertcat(CC,CS);
    A = zeros(size(C,1),3);
    for s=1:length(C)
        if C(s)==1
            A(s,:) = [0 0 .5]; % dark blue
        elseif C(s)==2
            A(s,:) = [.53 .83 .97]; % light blue
        elseif C(s)==3
            A(s,:) = [.7 .7 .7]; % grey
        elseif C(s)==4
            A(s,:) = [1 .60 .40]; % orange
        elseif C(s)==5
            A(s,:) = [.81 .07 .15]; % red
        elseif C(s)==6
            A(s,:) = [.99 .87 .09]; % yellow
        end
    end
    
end

if numROIs==34
h = scatter(xC,yC,S(pCort),A,'filled','MarkerEdgeColor',[.55 .55 .55],'LineWidth',0.7);
else
h = scatter(xC,yC,S(pCort),A(1:length(CC),:),'filled','MarkerEdgeColor',[.55 .55 .55],'LineWidth',0.7); hold on;
scatter(xS,yS,S(pSubcort),A(length(CC)+1:end,:),'o','LineWidth',1.5);
end
set(gcf,'color','w');
xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
title(sprintf('%s',whatNormalisation)); 
set(gca,'FontSize', 18);
cd ../../..
end
