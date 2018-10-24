function figure4AC(data2use, afterQC, options, Regeneratedata)
% make a figure for probe correlation from scratch
%1. regenerate data for different probe selection methods
if afterQC
    options.signalThreshold = 0.5;
    fileName = 'QC';

else
    options.signalThreshold = -1; % no QC filtering
    fileName = 'noQC';
end

options.VARfilter = false;
options.VARscale = 'normal';
options.VARperc = 50;
numIterRAND = 100;
saveOutput = true;

switch data2use
    case 'generate data'

        % files will be saved with standart name noQC

        if Regeneratedata
            options.probeSelections = {'Mean'};
            S2_probes(options);

            options.probeSelections = {'maxIntensity'};
            S2_probes(options);

            options.probeSelections = {'LessNoise'};
            S2_probes(options);

            options.probeSelections = {'PC'};
            S2_probes(options);

            options.probeSelections = {'DS'};
            S2_probes(options)

            options.probeSelections = {'Variance'};
            S2_probes(options)

            options.probeSelections = {'CV'};
            S2_probes(options)

            options.probeSelections = {'maxCorrelation_intensity'};
            S2_probes(options)

            options.probeSelections = {'maxCorrelation_variance'};
            S2_probes(options)
        end
        % load the initial data
        % created with there options
        % options.ExcludeCBandBS =  true;
        % options.useCUSTprobes = true;
        % options.updateProbes = 'reannotator';

        cd ('data/genes/processedData')
        load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
        if afterQC

            options.signalThreshold = 0.5; %  QC filtering
            options.useCUSTprobes = true;

            signalLevel = sum(noiseall,2)./size(noiseall,2);
            indKeepProbes = find(signalLevel>=options.signalThreshold);

            % sduplicated values calculated only on those probes that pass QC
            % threshold;
            [v, ind] = unique(DataTableProbe.EntrezID{1}(indKeepProbes));
            duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}(indKeepProbes), 1), ind);
            entrezSelected = DataTableProbe.EntrezID{1}(indKeepProbes);
            duplicate_value = unique(entrezSelected(duplicate_ind));
            percentage = (length(duplicate_value)/length(unique(entrezSelected)))*100;

        else

            [v, ind] = unique(DataTableProbe.EntrezID{1});
            duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
            duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));
            percentage = (length(duplicate_value)/length(unique(DataTableProbe.EntrezID{1})))*100;

        end

        format compact
        percentage

        % Load probes, selected using different methods that were just generated
        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXMean%s.mat', fileName))
        probes{1} = probeInformation;
        expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXLessNoise%s.mat', fileName))
        probes{2} = probeInformation;
        expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXPC%s.mat', fileName))
        probes{3} = probeInformation;
        expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXDS%s.mat', fileName))
        probes{4} = probeInformation;
        expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXVariance%s.mat', fileName))
        probes{5} = probeInformation;
        expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXCV%s.mat', fileName))
        probes{6} = probeInformation;
        expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxIntensity%s.mat', fileName))
        probes{7} = probeInformation;
        expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxCorrelation_intensity%s.mat', fileName))
        probes{8} = probeInformation;
        expression{8} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxCorrelation_variance%s.mat', fileName))
        probes{9} = probeInformation;
        expression{9} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
        %select only genes that had more than one probe available

        % do pairwise correlations for all non-random selection methods
        [~, indFilter] = intersect(probes{1}.EntrezID, duplicate_value);
        % filter genes that had multiple probes
        for k=1:length(probes)

            expression{k} = expression{k}(:,indFilter);
            probes{k}.EntrezID = probes{k}.EntrezID(indFilter);
        end

        avCorr = zeros(9,9);
        stdCorr = zeros(9,9);

        for i=1:9


            for j=i+1:9

                expr1 = expression{i};
                expr2 = expression{j};

                correlation = zeros(size(expr1,2),1);

                for g=1:size(expr1,2)
                    correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
                end

                avCorr(i,j) = mean(correlation);
                stdCorr(i,j) = std(correlation);
            end

        end
        avCorrfull = avCorr+avCorr';
        avCorrfull(logical(eye(size(avCorrfull)))) = 1;
        stdCorrfull = stdCorr+stdCorr';

        % loop over 100 random probe selections
        options.probeSelections = {'Random'};
        options.saveOutput = false;
        if afterQC
            options.signalThreshold = 0.5; % no QC filtering
            fileName = 'QC';
        else
            options.signalThreshold = -1; % no QC filtering
            fileName = 'noQC';
        end

        avCorrRAND = zeros(9,numIterRAND);
        stdCorrRAND = zeros(9,numIterRAND);

        %------------------------------------------------------------------------------
        % Load the data
        %------------------------------------------------------------------------------

        tic
        for iter=1:numIterRAND
            [expressionAll, probeInformation] = selectRANDprobes(options, DataTable, DataTableProbe,noiseall);
            probesRAND = probeInformation;
            expressionRAND = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
            % select only probes in other datasets
            expressionRAND = expressionRAND(:, indFilter);
            % for each iteration, correlate random probe celetion to other measures
            % and then take tha average after orver 100 runs


            for i=1:9

                expr1 = expression{i};
                expr2 = expressionRAND;

                correlation = zeros(size(expr1,2),1);

                for g=1:size(expr1,2)
                    correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
                end

                avCorrRAND(i,iter) = mean(correlation);
                stdCorrRAND(i,iter) = std(correlation);
            end
            fprintf('running iteration %d for all measures\n', iter)

        end
        toc
        % take the mean of 100 runs
        avCorrWrand = mean(avCorrRAND,2);
        stdCorrWrand = mean(stdCorrRAND,2);

        % add this row to the other matrix

        avCorrfullALL = [avCorrfull, avCorrWrand];
        avCorrWrand = [avCorrWrand;1];
        avCorrfullALL = [avCorrfullALL; avCorrWrand'];

        stdCorrfullALL = [stdCorrfull, stdCorrWrand];
        stdCorrWrand = [stdCorrWrand;0];
        stdCorrfullALL = [stdCorrfullALL; stdCorrWrand'];

        % reorder according to similarity
        R = BF_pdist(avCorrfullALL);
        [ord,R,keepers] = BF_ClusterReorder(avCorrfullALL,R);
        avCorrPlot = avCorrfullALL(ord, ord);
        stdCorrPlot = stdCorrfullALL(ord, ord);

        save(sprintf('probeCorrelationsRAND%s.mat',fileName),'avCorrPlot', 'stdCorrPlot', 'avCorrfullALL', 'avCorrRAND', 'stdCorrRAND');

        nice_cmap = [flipud(make_cmap('orangered',50,30,0))];

        figure; set(gcf,'Position',[300 300 1300 500])
        subplot(1,2,1) ; imagesc(avCorrPlot);
        set(gcf,'color','w');
        colormap(nice_cmap)
        caxis([0.5 1])
        colorbar
        tickNames = {'mean', 'signal proportion', 'PC', 'DS','variance', 'CV', 'max intensity', 'correlation intensity', 'correlation variance', 'random'};
        tickNamesORD = tickNames(ord);

        xticks([1 2 3 4 5 6 7 8 9 10])
        xticklabels(tickNamesORD);
        xtickangle(45)
        yticks([1 2 3 4 5 6 7 8 9 10])
        yticklabels(tickNamesORD);
        ytickangle(45)
        set(gca,'FontSize', 14)



        % make a plot with RNA seq
        options.probeSelections = {'RNAseq'};
        options.RNAseqThreshold = 0.2;
        options.useCUSTprobes = true;
        if Regeneratedata
            S2_probes(options)
        end

        % load data
        load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq%s.mat', fileName))
        probes{10} = probeInformation;
        expression{10} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
        % select genes in other versions and in this one  that overlap
        %entrezIDs = probes{1}.EntrezID(indFilter);

        [genesMultiple] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
        [genesMultipleRNAseq, indFilterRNA] = intersect(probes{10}.EntrezID, duplicate_value);

        % filter genes that had multiple probes

        expression{10} = expression{10}(:,indFilterRNA);

        [v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq);
        [v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple);
        % calculate correlation between each way of choosing a probe

        RNAcorr = zeros(9,1);
        stDEV = zeros(9,1);
        % calculate for RNAseq
        expr1 = expression{10}(:,indRNA);
        for j=1:9

            expr2 = expression{j}(:,indALL);
            correlation = zeros(size(expr1,2),1);

            for g=1:size(expr1,2)
                correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
            end


            RNAcorr(j) = mean(correlation);
            stDEV(j) = std(correlation);
        end

        % run a 100 runs of random selection and correlate each ro RNAseq
        options.probeSelections = {'Random'};
        options.saveOutput = false;
        if afterQC
            options.signalThreshold = 0.5; % no QC filtering
            fileName = 'QC';
        else
            options.signalThreshold = -1; % no QC filtering
            fileName = 'noQC';
        end

        %expr1 = expression{10}(:,indRNA); % this is RNAseq data
        avCorrRANDrna = zeros(numIterRAND,1);
        stdCorrRANDrna = zeros(numIterRAND,1);

        for iter = 1:numIterRAND
            [expressionAll, probeInformation] = selectRANDprobes(options, DataTable, DataTableProbe,noiseall);
            probesRAND = probeInformation;
            expressionRAND = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
            % select only probes in other datasets
            expressionRAND = expressionRAND(:, indFilter);% first filter as previous time
            % then filter based on RNAseq
            expr2 = expressionRAND(:,indALL);

            correlation = zeros(size(expr1,2),1);

            for g=1:size(expr1,2)
                correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
            end

            avCorrRANDrna(iter) = mean(correlation);
            stdCorrRANDrna(iter) = std(correlation);
            fprintf('running iteration %d for RNAseq\n', iter)
        end

        % get the mean for 100 iterations
        meanRNArand = mean(avCorrRANDrna);
        meanRNArandstd = mean(stdCorrRANDrna);

        RNAcorr = [RNAcorr; meanRNArand];
        stdevRNA = [stDEV; meanRNArandstd];

        [valS, indS] = sort(RNAcorr, 'descend');
        namesS = tickNames(indS);

        nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
        colors2use = nice_cmap([5 15 25 35 45 55 65 75 85 95],:);

        save(sprintf('probeCorrelationsRAND_RNAseq%s.mat', fileName), 'RNAcorr', 'stdevRNA', 'avCorrRANDrna', 'stdCorrRANDrna');

        subplot(1,2,2);
        set(gcf,'color','w');
        hold on
        for i = 1:length(valS)
            h=bar(i,valS(i));
            %errorbar(i,valS(i),stDEV(i),'.')
            if strcmp(namesS(i), 'CV')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'variance')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'random')
                set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'PC')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'mean')
                set(h,'FaceColor',[.95 .6 .6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'signal proportion')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'DS')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'max intensity')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'correlation intensity')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'correlation variance')

                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            end
        end
        hold off
        ylim([0.5 1]); ylabel('Spearman correlation'); %xlabel('Probe selection methods');
        %title(sprintf(' (%d)', length(indRNA)))
        xticks([1 2 3 4 5 6 7 8 9 10])
        %legend(sprintf('%d genes', length(indRNA)))
        xticklabels(namesS);
        xtickangle(45)
        set(gca,'FontSize', 14)

    case 'precomputed data'

        cd ('data/genes/processedData')
        load(sprintf('probeCorrelationsRAND%s.mat',fileName));
        nice_cmap = [flipud(make_cmap('orangered',50,30,0))];

        % reorder according to similarity
        R = BF_pdist(avCorrfullALL);
        [ord,R,keepers] = BF_ClusterReorder(avCorrfullALL,R);
        avCorrPlot = avCorrfullALL(ord, ord);

        figure; set(gcf,'Position',[300 300 1300 500])
        subplot(1,2,1); imagesc(avCorrPlot);
        set(gcf,'color','w');
        colormap(nice_cmap)
        caxis([0.5 1])
        colorbar
        tickNames = {'mean', 'signal proportion', 'PC', 'DS','variance', 'CV', 'max intensity', 'correlation intensity', 'correlation variance', 'random'};
        tickNamesORD = tickNames(ord);

        xticks([1 2 3 4 5 6 7 8 9 10])
        xticklabels(tickNamesORD);
        xtickangle(45)
        yticks([1 2 3 4 5 6 7 8 9 10])
        yticklabels(tickNamesORD);
        ytickangle(45)
        set(gca,'FontSize', 14)

        % RNAseq

        load(sprintf('probeCorrelationsRAND_RNAseq%s.mat', fileName));
        [valS, indS] = sort(RNAcorr, 'descend');
        namesS = tickNames(indS);

        subplot(1,2,2);
        set(gcf,'color','w');
        hold on
        for i = 1:length(valS)
            h=bar(i,valS(i));
            if strcmp(namesS(i), 'CV')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'variance')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'random')
                set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'PC')
                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'mean')
                set(h,'FaceColor',[.95 .6 .6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'signal proportion')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'DS')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'max intensity')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'correlation intensity')
                set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            elseif strcmp(namesS(i), 'correlation variance')

                set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
            end
        end
        hold off
        ylim([0.5 1]); ylabel('Spearman correlation');
        xticks([1 2 3 4 5 6 7 8 9 10])
        xticklabels(namesS);
        xtickangle(45)
        set(gca,'FontSize', 14)
end
cd ../../..
end
