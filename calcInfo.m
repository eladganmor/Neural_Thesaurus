%calculate information where data is numOfReps x numOfStim matrix, and
%entries are taken from {1...range};

%   Copyright 2015 Elad Ganmor
function [info h hCond] = calcInfo(data)
    range = max(data(:));
    [numOfReps, numOfStim] = size(data);
    condProb = zeros(range,numOfStim);
    for s=1:numOfStim
        counts = hist(data(:,s),1:range);
        condProb(:,s) = counts/numOfReps;
    end
    condProb = condProb + eps;
    hCond = mean(diag(-condProb'*log2(condProb)));%conditional entropy
    prob = sum(condProb,2)/numOfStim;
    h = -prob'*log2(prob);%entropy
    info = h - hCond;
end
    
