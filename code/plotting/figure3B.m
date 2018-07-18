function figure3B()
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

signalThreshold = [-1 0.5];
figure; colors = [.96 .63 .55; 1 .46 .22];
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
        %howMany = length(A);
        if isMore
            if length(A)>numProbes
                m=m+1;
                %PL{p} = A;
                r = NaN(length(A));
                for k=1:length(A)
                    for l=k+1:length(A)
                        r(k,l) = corr(Expressionall2(A(k),:)', Expressionall2(A(l),:)', 'type', 'Spearman');
                    end
                    %IND2(w) = A(k);
                    w=w+1;
                end
                t=r(:);
                t(isnan(t)) = [];
                corVal(p) = mean(t);
                %w=w+length(A);
            end
        else
            if length(A)==numProbes
                m=m+1;
                %PL{p} = A;
                r = NaN(length(A));
                for k=1:length(A)
                    for l=k+1:length(A)
                        r(k,l) = corr(Expressionall2(A(k),:)', Expressionall2(A(l),:)', 'type', 'Spearman');
                    end
                    %IND2(w) = A(k);
                    w=w+1;
                end
                t=r(:);
                t(isnan(t)) = [];
                corVal(p) = mean(t);
                %w=w+length(A);
            end
            
            
        end
    end
    
    
    inds = ~isnan(corVal);
    multind = find(inds==1);
    corMult{i} = corVal(multind);
    
    perc = length(find(corMult{i}<0.3))/length(corMult{i});
    
    histogram(corMult{i}, 50,'EdgeColor',[.6 .6 .6],...
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
cd ../../..
end
