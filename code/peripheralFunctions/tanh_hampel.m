function [x_hat,locEstimate,scaleEstimate] = tanh_hampel(x,a,b,c)

% Do all normalization relative to the median of the distribution
goodR = ~isnan(x);
xGood = x(goodR);
% xRel = x(goodR) - median(x(goodR));

if nargin < 2
    % set 70% of data within a of median:
    % (as https://books.google.com.au/books?id=1Wpx25D8qOwC&pg=PA328&lpg=PA328&dq=hampel+tanh&source=bl&ots=9xO_7Snrb6&sig=bINDZsHbAZq6nK5FKBE9Lx6HORo&hl=en&sa=X&ved=0CB0Q6AEwAGoVChMI5dDpw5SjyAIVhSWUCh3Vqw4M#v=onepage&q=hampel%20tanh&f=false)
    a = withinXofmedian(0.7);
    b = withinXofmedian(0.85);
    c = withinXofmedian(0.95);
end

[locEstimate,fval] = fminsearch(@sumInfluence,median(xGood));
scaleEstimate = sqrt(mean(getInfluence(locEstimate).^2));

%-------------------------------------------------------------------------------
% locEstimate = mean(infl);
% scaleEstimate = std(infl);
%-------------------------------------------------------------------------------
x_hat = zeros(size(x));
x_hat(goodR) = 0.5*(tanh(0.1*((xGood-locEstimate)/scaleEstimate))+1);
x_hat(~goodR) = NaN;
% x_hat = UnityRescale(x_hat);

% if corr(x,x_hat)<0.99
    % subplot(3,1,1)
    % histogram(x,25);
    % subplot(3,1,2)
    % histogram(x_hat,25);
    % subplot(3,1,3);
    % plot(x,x_hat,'.k')
    % keyboard
% end

%-------------------------------------------------------------------------------
function f = withinXofmedian(prop)
    f = quantile(abs(xGood-median(xGood)),prop);
end
%-------------------------------------------------------------------------------
function Suminfl = sumInfluence(xLoc)
    infl = getInfluence(xLoc);
    % Sum of influence function:
    Suminfl = abs(sum(infl));
    % plot(xRel,infl,'.k')
end
%-------------------------------------------------------------------------------
function infl = getInfluence(xLoc)
    xRel = xGood - xLoc;

    %-------------------------------------------------------------------------------
    % Determine the influence function:
    infl = zeros(length(xRel),1);

    % Points close to zero (within a) count the most:
    r1 = abs(xRel)<a;
    infl(r1) = xRel(r1);
    % Points between a and b count as a:
    r2 = abs(xRel)>=a & abs(xRel)<b;
    infl(r2) = a*sign(xRel(r2));
    % Points between b and c count as linear decreasing amount:
    r3 = abs(xRel)>=b & abs(xRel)<c;
    infl(r3) = a*sign(xRel(r3)).*(c-abs(xRel(r3)))/(c-b);
    % Points greater than c count as zero:
    % --
end
%-------------------------------------------------------------------------------
function xhat = UnityRescale(x)
    % Linearly rescale a data vector to unit interval:
    goodVals = ~isnan(x); % only work with non-NaN data
    if ~any(goodVals) % all NaNs:
        xhat = x; return
    end
    minX = min(x(goodVals));
    maxX = max(x(goodVals));
    if minX==maxX
        % There is no variation -- set all to zero -- attempts to rescale will blow up
        xhat = x;
        xhat(goodVals) = 0;
    else
        % There is some variation -- rescale to unit interval
        xhat = (x-minX)/(maxX-minX);
    end
end

end
