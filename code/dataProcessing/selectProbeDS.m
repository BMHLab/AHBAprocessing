% Dscalculation based on overlapping structures for probe selection
function [indProbe,avCorr] = selectProbeDS(EntrezID, DataTable, indKeepProbes)

% for each pair of subjects find overlapping samples, average data and
% calculate correlation between each probe for the same gene
ugenes = unique(EntrezID);

probeS = cell(6,6); %,length(ugenes));
mCorrs = cell(6,6);

for s1=1:6
    samp1 = DataTable.SampleID{s1};
    usamp1 = unique(samp1);
    % average expression for the same region
    for r1=1:length(usamp1)
        indr1 = samp1==usamp1(r1);
        expS1(:,r1) = mean(DataTable.Expression{s1}(indKeepProbes,indr1),2);
    end
    
    for s2=s1+1:6
        samp2 = DataTable.SampleID{s2};
        usamp2 = unique(samp2);
        % average expression for the same region
        for r2=1:length(usamp2)
            indr2 = samp2==usamp2(r2);
            expS2(:,r2) = mean(DataTable.Expression{s2}(indKeepProbes,indr2),2);
        end
        
        % find intersect between structures for S1 and S2
        [~, rindS1, rindS2] = intersect(usamp1, usamp2);
        
        Pr =  zeros(length(ugenes),1);
        maxCor = zeros(length(ugenes),1);
        
        for g = 1:length(ugenes)
            % find unique genes
            gprobes = find(EntrezID==ugenes(g));
            if length(gprobes)>1
                C = zeros(length(gprobes),1);
                for d = 1:length(gprobes)
                    
                    C(d) = corr(expS1(gprobes(d),rindS1)', expS2(gprobes(d),rindS2)', 'type', 'Spearman');
                    
                end
                [valM,indM] = max(C);
                Pr(g) = gprobes(indM);
               % Pr(g) = indM;
                maxCor(g) = valM;
            else
                Pr(g) = gprobes;
                maxCor(g) = corr(expS1(gprobes,rindS1)', expS2(gprobes,rindS2)', 'type', 'Spearman');
            end
            
            
        end
        probeS{s1,s2} = Pr;
        mCorrs{s1,s2} = maxCor;
        
    end

end

P = horzcat(probeS{1,2}, probeS{1,3}, probeS{1,4}, probeS{1,5}, probeS{1,6}, ...
    probeS{2,3}, probeS{2,4}, probeS{2,5}, probeS{2,6}, probeS{3,4}, probeS{3,5}, ...
    probeS{3,6}, probeS{4,5}, probeS{4,6}, probeS{5,6});

C = horzcat(mCorrs{1,2}, mCorrs{1,3}, mCorrs{1,4}, mCorrs{1,5}, mCorrs{1,6}, ...
    mCorrs{2,3}, mCorrs{2,4}, mCorrs{2,5}, mCorrs{2,6}, mCorrs{3,4}, mCorrs{3,5}, ...
    mCorrs{3,6}, mCorrs{4,5}, mCorrs{4,6}, mCorrs{5,6});


indProbe = mode(P,2);
avCorr = mean(C,2);
end

