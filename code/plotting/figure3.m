function figure3()

cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
%% 3A
fprintf('Making figure 3A\n'); 
signalThresholds = 0:0.01:1;
% count how many genes/probes are there and how many are excluded with
% signal threshold of 0.5
probes = DataTableProbe.ProbeName{1}; 
genes = DataTableProbe.EntrezID{1}; 

[uGenes, uGenesIND] = unique(genes);
fprintf('%d unique probes\n', length(probes))
fprintf('%d unique genes\n', length(uGenes))

%% find repeating entrezIDs and calculate variances for them, then take a probe of max variance.

% % ------------------------------------------------------------------------------
% First, find probes that have very noisy data and remove them from consideration
% Threshold for removing those proges is defined as the percentage of
% samples a probe has expression higher that background
% % ------------------------------------------------------------------------------
nrProbes = zeros(length(signalThresholds),1); 
nrGenes = zeros(length(signalThresholds),1); 

for j=1:length(signalThresholds)
    
signalThreshold = signalThresholds(j);     
    
signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=signalThreshold);
indRemoveProbes = find(signalLevel<signalThreshold);
keepProbes = signalLevel>=signalThreshold;

% remove selected probes from data and perform other calculations only on
% non-noisy probes
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);

nrProbes(j) = length(ProbeName); 
nrGenes(j) = length(unique(EntrezID)); 
end

fig=figure;
set(fig,'Position', [100, 100, 1500, 400]);
set(fig,'defaultAxesColorOrder',[.92 .35 .24; 0 .5 .5]);
subplot(1,3,1);

yyaxis left
plot(signalThresholds, nrGenes, '-o', ...
    'Color', [.6 .6 .6], 'LineWidth',1,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[.6 .6 .6],...
    'MarkerFaceColor',[.96 .63 .55]);
xlabel('Intensity-based filtering threshold'); 
ylabel('Number of genes')
yticklabels([0 5000 10000 15000 20000 25000])
yticks([0 5000 10000 15000 20000 25000])
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
set(gcf,'color','w');
hold on; 

yyaxis right
plot(signalThresholds, nrProbes, '-o', ...
    'Color', [.6 .6 .6], 'LineWidth',1,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[.6 .6 .6],...
    'MarkerFaceColor',[.59 .87 .82]);
ylabel('Number of probes','FontSize', 14)
yticklabels([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
yticks([0 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000])
%legend({'Number of genes', 'Number of probes'}); 
set(gca,'FontSize', 16)
box off
h = legend({'Number of genes', 'Number of probes'}); 
set(h,'fontsize',20)

%% 3B
fprintf ('Making figure 3B\n')
signalThreshold = [-1 0.5];
colors = [.96 .63 .55; 1 .46 .22];
corMult = cell(length(signalThreshold),1);
numProbes = 1;
isMore = true;
for i=1:length(signalThreshold)
    
    signalLevel = sum(noiseall,2)./size(noiseall,2);
    indKeepProbes = find(signalLevel>=signalThreshold(i));
    
    
    [v, ind] = unique(DataTableProbe.EntrezID{1}(indKeepProbes));
    entrezID = DataTableProbe.EntrezID{1}(indKeepProbes);
    Expressionall2 = Expressionall((indKeepProbes),:);
    
    m=0; w=1;
    corVal = nan(length(ind),1);
    for p=1:length(ind)
        A = find(entrezID==v(p));
        if isMore
            if length(A)>numProbes
                m=m+1;
                r = NaN(length(A));
                for k=1:length(A)
                    for l=k+1:length(A)
                        r(k,l) = corr(Expressionall2(A(k),:)', Expressionall2(A(l),:)', 'type', 'Spearman');
                    end
                    w=w+1;
                end
                t=r(:);
                t(isnan(t)) = [];
                corVal(p) = mean(t);
            end
        else
            if length(A)==numProbes
                m=m+1;
                r = NaN(length(A));
                for k=1:length(A)
                    for l=k+1:length(A)
                        r(k,l) = corr(Expressionall2(A(k),:)', Expressionall2(A(l),:)', 'type', 'Spearman');
                    end
                    w=w+1;
                end
                t=r(:);
                t(isnan(t)) = [];
                corVal(p) = mean(t);
            end
            
            
        end
    end
    
    
    inds = ~isnan(corVal);
    multind = find(inds==1);
    corMult{i} = corVal(multind);
    
    perc = length(find(corMult{i}<0.3))/length(corMult{i});
    
    subplot(1,3,2); histogram(corMult{i}, 50,'EdgeColor',[.6 .6 .6],...
        'FaceColor',colors(i,:));
    
    set(gcf,'color','w'); hold on;
end

xlabel('Average correlation between probes')
ylabel('Number of genes')
set(gca,'FontSize', 16)
legendText{1} = sprintf('Before intensity-based filtering (%d genes)', length(corMult{1}));
legendText{2} = sprintf('After intensity-based filtering (%d genes)', length(corMult{2}));
h = legend(legendText{1},legendText{2});
set(h,'fontsize',20)
xticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1])

box off

%% 3C
fprintf ('Making figure 3C\n')
clearvars -except DataTableProbe Expressionall noiseall 

for j=1:2
    % after filtering
    if j==1
        doFilter = true;
    elseif j==2
        doFilter = false;
    end
    
    
    signalLevel = sum(noiseall,2)./size(noiseall,2);
    indKeepProbes = find(signalLevel>=0.5);
    
    EntrezIDfiltr = DataTableProbe.EntrezID{1,1}(indKeepProbes);
    Expressionallfiltr = Expressionall(indKeepProbes,:);
    
    
    [v, ind] = unique(EntrezIDfiltr);
    
    duplicate_ind = setdiff(1:size(EntrezIDfiltr), ind);
    duplicate_value = unique(EntrezIDfiltr(duplicate_ind));
    
    
    if doFilter
        Ent = EntrezIDfiltr;
        Expr = Expressionallfiltr;
    else
        Ent = DataTableProbe.EntrezID{1,1};
        Expr = Expressionall;
    end
    
    
    for i=1:length(duplicate_value)
        
        
        A = find(Ent==duplicate_value(i));
        
        if length(A)>1
            
            rp = NaN(length(A));
            for k=1:length(A)
                for l=k+1:length(A)
                    rp(k,l) = corr(Expr(A(k),:)', Expr(A(l),:)', 'type', 'Spearman');
                end
            end
            t=rp(:);
            t(isnan(t)) = [];
  
            
            corVal(i,j) = mean(t);
            
        end
        
    end
end

[pall,hall,statsall] = ranksum(corVal(:,1),corVal(:,2));

% select only genes that were affected by this filtering - the ones where
% correlation between probes changed.
for i=1:size(corVal,1)
    if corVal(i,1) == corVal(i,2)
        indRem(i) = i;
    end
end

indRem(indRem==0) = [];
C = corVal;
C(indRem,:) = [];

mafter = mean(C(:,1));
mbefore = mean(C(:,2));

subplot(1,3,3);
histogram(C(:,2), 25,'EdgeColor',[.6 .6 .6],...
    'FaceColor',[.96 .63 .55]);
xlabel('Average correlation between probes')
ylabel('Number of genes')
set(gcf,'color','w'); hold on;

histogram(C(:,1), 30,'EdgeColor',[.6 .6 .6],...
    'FaceColor',[1 .46 .22]);

legendText{2} = sprintf('After intensity-based filtering (%d genes)', size(C,1));
legendText{1} = sprintf('Before intensity-based filtering (%d genes)', size(C,1));
set(gca,'FontSize', 16)
h = legend(legendText{1},legendText{2});
set(h,'fontsize',20)
xticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1])

box off

[psubset,hsubset,statssubset] = ranksum(C(:,1),C(:,2));
cd ../../..
end

