function [expPlot, parcelCoexpression, correctedCoexpression, Residuals, distExpVect, c, parcelExpression] = calculateCoexpression(distRegions, selectedGenes, DSvalues, ROIs,nROIs, Fit, correctDistance, xrange, doPlotCGE, doPlotResiduals, ROIind, how2mean, percentDS)

if nargin<9
    doPlotCGE = false;
    doPlotResiduals = false;
end


if nargin<10
    doPlotResiduals = false;
end

if length(ROIind)==6
    subjIND = vertcat(ROIind{1}(:,2), ROIind{2}(:,2), ROIind{3}(:,2), ROIind{4}(:,2), ROIind{5}(:,2), ROIind{6}(:,2));
elseif length(ROIind)==2
    subjIND = vertcat(ROIind{1}(:,2), ROIind{2}(:,2));
end

selectedGenesN = selectedGenes;
parcelExpression = zeros(length(nROIs), size(selectedGenesN,2));

if strcmp(how2mean, 'meanSamples')

    for samp=1:length(nROIs)
        
        A = ROIs == nROIs(samp);
        if length(find(A))>1
            parcelExpression(samp,:) = nanmean(selectedGenesN(A==1,:));
        elseif length(find(A))==1
            parcelExpression(samp,:) = selectedGenesN(A==1,:);
        else
            parcelExpression(samp,:) = NaN;
        end
        
    end
    
elseif strcmp(how2mean, 'meanSubjects')
    parcelExpressionSubj = zeros(6, length(nROIs), size(selectedGenesN,2));
    
    for samp=1:length(nROIs)
        for subj=1:size(ROIind,2)
            isSample = subjIND==subj;
            isSubject = ROIs == nROIs(samp);
            
            A = zeros(length(subjIND),1);
            A(isSample==1 & isSubject==1) = 1;
            
            if length(find(A))>1
                parcelExpressionSubj(subj,samp,:) = nanmean(selectedGenesN(A==1,:));
            elseif length(find(A))==1
                parcelExpressionSubj(subj,samp,:) = selectedGenesN(A==1,:);
            else
                parcelExpressionSubj(subj,samp,:) = NaN;
            end
            
        end
    end
    parcelExpression = squeeze(nanmean(parcelExpressionSubj,1));
end

% select DS genes for CGE calculation - all other data is the same
if percentDS<100
selectedGenesParc = parcelExpression(:,DSvalues(:,1)); % take genes with highest DS values
else
selectedGenesParc = parcelExpression;
end
parcelCoexpression = corr(selectedGenesParc', 'type', 'Spearman'); % calculate sample-sample coexpression
parcelCoexpression(logical(eye(size(parcelCoexpression)))) = NaN; % replace diagonal with NaN

distExpVect(:,1) = distRegions(:);
distExpVect(:,2) = parcelCoexpression(:);

Dvect = distExpVect(:,1);
Rvect = distExpVect(:,2);

%----------------------------------------------------------------------------------
% Fit distance correction according to a defined rule
%----------------------------------------------------------------------------------

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

numSamples = size(distRegions,1);
correctedCoexpression = reshape(Residuals,[numSamples, numSamples]);

%----------------------------------------------------------------------------------
% Plot corrected ROI-ROI coexpression matrix
%----------------------------------------------------------------------------------

if correctDistance
    expPlot = correctedCoexpression;
else
    expPlot = parcelCoexpression;
end

end
