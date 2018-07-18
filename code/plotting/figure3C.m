function figure3C()
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
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
            
            r = NaN(length(A));
            for k=1:length(A)
                for l=k+1:length(A)
                    r(k,l) = corr(Expr(A(k),:)', Expr(A(l),:)', 'type', 'Spearman');
                end
                
                
            end
            t=r(:);
            t(isnan(t)) = [];
            
            %corVal(i,j) = median(t);
            
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

figure;
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
