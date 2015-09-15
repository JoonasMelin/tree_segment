function [rangeMin, rangeMax] = getValueRange(img, includeRange)
[numVals, binMean] = hist(medfilt1(img(:),10), 20000);
cumVals = medfilt1(cumsum(numVals));
minVal = min(cumVals);
maxVal = max(cumVals);
rangeVal = maxVal-minVal;
minLim = minVal+(rangeVal*(1-includeRange)/2);
maxLim = maxVal-(rangeVal*(1-includeRange)/2);

minInd = find(cumVals >= minLim, 1, 'first');
maxInd = find(cumVals >= maxLim, 1, 'first');



if minInd > 0
    rangeMin = binMean(minInd-1);
else
    rangeMin = binMean(minInd);
end

if maxInd < length(binMean)
    rangeMax = binMean(maxInd+1);
else
    rangeMax = binMean(maxInd);
end

% figure(1);
% plot(cumVals);
% hold on;
% plot(minInd, cumVals(minInd), 'xb', 'LineWidth', 3);
% plot(maxInd, cumVals(maxInd), 'xr', 'LineWidth', 3);
% hold off;