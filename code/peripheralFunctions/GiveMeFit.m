function [f_handle,Stats,c] = GiveMeFit(xData,yData,fitType,suppressFigure)
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-11-20
% ------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (~isnan(xData) & ~isnan(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

if nargin < 4
    suppressFigure = 0;
end

switch fitType
case 'decayEta0'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[200,1]);
    f = fittype('A*x.^-n','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to decayEta0 failed')
    end
    f_handle = @(x) c.A.*x.^(-c.n);

case 'decayEta'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[2,1,-1]);
    f = fittype('A*x.^-n + B','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to decayEta failed')
    end
    f_handle = @(x) c.A.*x.^(-c.n) + c.B;

case 'exp'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,0]);
    f = fittype('A*exp(-n*x) + B','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to exp failed')
    end
    f_handle = @(x) c.A.*exp(-c.n*x) + c.B;

case 'exp0'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1]);
    f = fittype('A*exp(-n*x)','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to exp failed')
    end
    f_handle = @(x) c.A.*exp(-c.n*x);

case 'exp1'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0]);
    f = fittype('exp(-n*x) + B','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to exp1 failed')
    end
    f_handle = @(x) exp(-c.n*x) + c.B;

case 'exp_1_0'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',1);
    f = fittype('exp(-n*x)','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to exp_1_0 failed')
    end
    f_handle = @(x) exp(-c.n*x);

    if ~suppressFigure
        figure('color','w'); box('on');
        plot(xData,yData,'.k'); hold on; plot(xData,f_handle(xData),'.r');
        title(sprintf('Fit exponential to %u points with eta = %.4g\n',length(xData),c.n),'interpreter','none')
    end
    fprintf(1,'Fit exponential with eta = %.4g\n',c.n);

case 'decay0'
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',200);
    f = fittype('A/x','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to decay0 failed')
    end
    f_handle = @(x) c.A./x;

case 'decay' % 1/x decay

    % f = fittype('A*exp(B*x) + C','options',s);
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[500, min(yData)]);
    f = fittype('A/x + B','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
    catch
        error('Fit to decay failed')
    end
    f_handle = @(x) c.A./x + c.B;

case 'linear'
    [c,Stats] = fit(xData,yData,'poly1');
    % [p, Stats] = polyfit(xData,yData,1);
    Gradient = c.p1; Intercept = c.p2;
    f_handle = @(x) Gradient*x + Intercept;

case 'robust_linear'
    [p, Stats] = robustfit(xData,yData);
    Gradient = p(2); Intercept = p(1);
    f_handle = @(x) Gradient*x + Intercept;

otherwise
    error('Unknown fit type: ''%s''',fitType);
end

end
