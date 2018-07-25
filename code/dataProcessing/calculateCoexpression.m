function [expPlot, parcelCoexpression, correctedCoexpression, Residuals, distExpVect, distPlot, c, parcelExpression] = calculateCoexpression(sampleDistances, selectedGenes, DSvalues, W, ROIs,nROIs, Fit, correctDistance, resolution, xrange, doPlotCGE, doPlotResiduals, ROIind, how2mean)

if nargin<12
    doPlotCGE = false;
    doPlotResiduals = false;
end


if nargin<13
    doPlotResiduals = false;
end

subjIND = vertcat(ROIind{1}(:,2), ROIind{2}(:,2), ROIind{3}(:,2), ROIind{4}(:,2), ROIind{5}(:,2), ROIind{6}(:,2)); 
selectedGenesN = selectedGenes(:,DSvalues(:,1)); % take genes with highest DS values
switch resolution
    case 'sample'
        
        sampleCoexpression = corr(selectedGenesN', 'type', 'Spearman'); % calculate sample-sample coexpression
        
        
        sampleCoexpression(logical(eye(size(sampleCoexpression)))) = NaN; % replace diagonal with NaN
        distExpVect(:,1) = sampleDistances(:); % make a vector for distances
        distExpVect(:,2) = sampleCoexpression(:); % make a vector for coexpression values
        distExpVect(any(isnan(distExpVect), 2), :) = [];
        
        Dvect = distExpVect(:,1);
        Rvect = distExpVect(:,2);
        
    case 'ROI'
        parcelExpression = zeros(length(W), size(selectedGenesN,2));
        
        if strcmp(how2mean, 'meanSamples')
            
            
            for samp=1:length(W)
                
                A = ROIs == W(samp);
                if length(find(A))>1
                    parcelExpression(samp,:) = nanmean(selectedGenesN(A==1,:));
                else
                    parcelExpression(samp,:) = selectedGenesN(A==1,:);
                end
                
            end
            
        elseif strcmp(how2mean, 'meanSubjects')
            parcelExpressionSubj = zeros(6, length(W), size(selectedGenesN,2));
            
            for samp=1:length(W)
                for subj=1:size(ROIind,2)
                    isSample = subjIND==subj; 
                    isSubject = ROIs == W(samp);
                    
                    A = zeros(length(subjIND),1); 
                    A(isSample==1 & isSubject==1) = 1; 

                    if length(find(A))>1
                        parcelExpressionSubj(subj,samp,:) = nanmean(selectedGenesN(A==1,:));
                    elseif isempty(find(A, 1))
                        parcelExpressionSubj(subj,samp,:) = NaN;
                    else
                        parcelExpressionSubj(subj,samp,:) = selectedGenesN(A==1,:);
                    end

                end
            end
            parcelExpression = squeeze(nanmean(parcelExpressionSubj,1)); 

        end
        
        [sROIs, ind] = sort(ROIs);
        distancesSorted = sampleDistances(ind, ind);
        parcelDistances = zeros(length(W),length(W));
        
        for sub=1:length(W)
            for j=1:length(W)
                
                A = sROIs == W(sub);
                isSubject = sROIs == W(j);
                
                D = distancesSorted(A,isSubject);
                parcelDistances(sub,j) = nanmean(D(:));
                
            end
        end
        
        p = setdiff(nROIs, ROIs);
        parcelCoexpression = corr(parcelExpression', 'type', 'Spearman'); % calculate sample-sample coexpression
        parcelCoexpression(logical(eye(size(parcelCoexpression)))) = NaN; % replace diagonal with NaN
        
        if ~isempty(p)
            
            parcelDistances = insertrows(parcelDistances,NaN,p);
            parcelDistances = insertrows(parcelDistances.', NaN,p).' ; % insert columns
            
            parcelCoexpression = insertrows(parcelCoexpression,NaN,p);
            parcelCoexpression = insertrows(parcelCoexpression.', NaN,p).' ; % insert columns
            
        end
        
        distExpVect(:,1) = parcelDistances(:);
        distExpVect(:,2) = parcelCoexpression(:);
        
        
        Dvect = distExpVect(:,1);
        Rvect = distExpVect(:,2);
        %distExpVect(any(isnan(distExpVect), 2), :) = [];
end
%----------------------------------------------------------------------------------
% Fit distance correction according to a defined rule
%----------------------------------------------------------------------------------
%

% select values with distance <100
if strcmp(Fit{1}, 'linear') || strcmp(Fit{1}, 'exp') || strcmp(Fit{1}, 'exp_1_0') || strcmp(Fit{1}, 'decay') || strcmp(Fit{1}, 'exp0') || strcmp(Fit{1}, 'exp1')
    
    [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
else
    c=0;
end

% plot original coexpression-distance .
if doPlotCGE
    [xThresholds,yMeans] = BF_PlotQuantiles(distExpVect(:,1),distExpVect(:,2),26,1,1);
    xlabel('Distance between regions (mm)'); ylabel('Correlated gene expression');
    ylim([-1 1]); xlim(xrange);
    set(gca,'fontsize',18)
end
switch Fit{1}
    
    case 'linear'
        FitCurve = c.p1*Dvect + c.p2;
    case 'exp'
        FitCurve = c.A*exp(-c.n*Dvect) + c.B;
    case 'exp_1_0'
        FitCurve = exp(-c.n*Dvect);
    case 'decay'
        FitCurve = c.A/Dvect + c.B;
        Residuals = Rvect' - FitCurve;
    case 'exp0'
        FitCurve = c.A.*exp(-c.n*Dvect);
    case 'exp1'
        FitCurve = exp(-c.n*Dvect) + c.B;
    otherwise
        Y = discretize(distExpVect(:,1),xThresholds);
        Residuals = zeros(length(Y),1);
        for val=1:length(Y)
            if ~isnan(distExpVect(val,2))
                Residuals(val) = distExpVect(val,2) - yMeans(Y(val));
            else
                Residuals(val) = NaN;
            end
            
            
        end
        
end

if strcmp(Fit{1}, 'linear') || strcmp(Fit{1}, 'exp') || strcmp(Fit{1}, 'exp_1_0') || strcmp(Fit{1}, 'decay') || strcmp(Fit{1}, 'exp0') || strcmp(Fit{1}, 'exp1')
    if doPlotCGE
        hold on; scatter(distExpVect(:,1),FitCurve,3, '.', 'r');
        xlim(xrange)
    end
    % get residuals
    Residuals = Rvect - FitCurve;
    if doPlotResiduals
        BF_PlotQuantiles(distExpVect(:,1),nonzeros(Residuals(:)),26,1,1);   xlabel('Distance between regions (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',18)
        xlim(xrange)
    end
else
    if doPlotResiduals
        BF_PlotQuantiles(distExpVect(:,1),Residuals(:),26,1,1);   xlabel('Distance between regions (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',18)
        xlim(xrange)
    end
end

%----------------------------------------------------------------------------------
% Plot corrected sample - sample coexpression matrix
%----------------------------------------------------------------------------------
switch resolution
    case 'sample'
        numSamples = size(sampleDistances,1);
        % add NaNs to diagonal for reshaping
        Idx=linspace(1, size(sampleDistances,1)*size(sampleDistances,1),size(sampleDistances,1));
        c=false(1,length(Residuals)+length(Idx));
        c(Idx)=true;
        nResiduals = nan(size(c));
        nResiduals(~c) = Residuals;
        nResiduals(c) = NaN;
        
        correctedCoexpression = reshape(nResiduals,[numSamples, numSamples]);
    case 'ROI'
        numSamples = size(parcelDistances,1);
        sampleDistances = parcelDistances;
        correctedCoexpression = reshape(Residuals,[numSamples, numSamples]);
        
end

% if doPlot
%     figure; imagesc(correctedCoexpression); caxis([-1,1]);title('Corrected coexpression');
%     colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
% end

%----------------------------------------------------------------------------------
% Plot corrected ROI-ROI coexpression matrix
%----------------------------------------------------------------------------------
switch resolution
    case 'sample'
        
        parcelCoexpression = zeros(length(W),length(W));
        parcelDistances = zeros(length(W),length(W));
        
        [sROIs, ind] = sort(ROIs);
        correctedCoexpressionSorted = correctedCoexpression(ind, ind);
        
        sampleDistances(logical(eye(size(sampleDistances)))) = NaN;
        distancesSorted = sampleDistances(ind, ind);
        
        
        sampleCoexpression(logical(eye(size(sampleCoexpression)))) = NaN;
        coexpressionSorted = sampleCoexpression(ind, ind);
        %         if doPlot
        %             figure; subplot(1,25,[1 2]); imagesc(sROIs);
        %             subplot(1,25,[3 25]); imagesc(coexpressionSorted); caxis([-1,1]); title('Corrected coexpression sorted samples');
        %             colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
        %             set(gca,'xtick',[])
        %             set(gca,'xticklabel',[])
        %             set(gca,'ytick',[])
        %             set(gca,'yticklabel',[])
        %         end
        
        isNormP = zeros(length(W),length(W));
        isNormH = zeros(length(W),length(W));
        
        for sub=1:length(W)
            for j=1:length(W)
                
                A = sROIs == W(sub);
                isSubject = sROIs == W(j);
                %for corrected
                if correctDistance == true
                    
                    P = correctedCoexpressionSorted(A, isSubject);
                    
                else
                    
                    P = coexpressionSorted(A, isSubject);
                    
                end
                D = distancesSorted(A,isSubject);
                parcelDistances(sub,j) = mean(D(:));
                parcelCoexpression(sub,j) = mean(P(:));
                
            end
        end
        
end

if correctDistance == true && strcmp(resolution, 'ROI')
    expPlot = correctedCoexpression;
else
    expPlot = parcelCoexpression;
end

distPlot = parcelDistances;

if strcmp(resolution, 'sample')
    p = setdiff(nROIs, ROIs);
    if ~isempty(p)
        
        expPlot = insertrows(expPlot,NaN,p);
        expPlot = insertrows(expPlot.', NaN,p).' ; % insert columns
        
        distPlot = insertrows(distPlot,NaN,p);
        distPlot = insertrows(distPlot.', NaN,p).' ; % insert columns
        
    end
end
% if doPlot
%     figure('color','w'); box('off');
%     imagesc(expPlot); caxis([-1,1]);
%     if correctDistance
%         title('Corrected parcellation coexpression ROIs');
%     else
%         title('Non - corrected parcellation coexpression ROIs');
%     end
%     colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
% end
end
