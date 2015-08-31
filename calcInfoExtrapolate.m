%calculate information where data is numOfReps x numOfStim matrix, and
%extrapolate to infinite data (Strong et al. PRL 1998);
%Input: 1. Array of N_repeats x N_stim of responses (as decimal numbers)
%       2. Number of repeats to perform for the extrapolation (more repeats better accuracy, but slower)
%Output: extrapolated information, along with extrapolated entropy and
%conditional entropy

%   Copyright 2015 Elad Ganmor
function [infoExtrapolate hExtrapolate hCondExtrapolate] = calcInfoExtrapolate(data,repeats)
    if nargin<2
        repeats = 1;
    end
    dataFract = 0.5:0.05:0.9;%for extrapolation 
    [numOfReps, numOfStim] = size(data);
    info = zeros(repeats,length(dataFract)); %mutual information
    h = zeros(repeats,length(dataFract)); %entropy
    hGiven = zeros(repeats,length(dataFract)); %conditional entropy
    
    parfor i=1:repeats
        %use tmp arrays only necessary for parfor
        tmpI = zeros(1,length(dataFract));
        tmpH = zeros(1,length(dataFract));
        tmpHGiven = zeros(1,length(dataFract));
        
        for dataFractInd = 1:length(dataFract)
            p = randperm(numOfReps); %choose subset of repeats
            chosenReps = p(1:round(dataFract(dataFractInd)*numOfReps));
            [tmpI(dataFractInd) tmpH(dataFractInd) tmpHGiven(dataFractInd)] = calcInfo(data(chosenReps,:)); %calculate information
        end
        
        %results for this repeat
        info(i,:) = tmpI;
        h(i,:) = tmpH;
        hGiven(i,:) = tmpHGiven;
    end
    
    %use full data set (no point in repeating this calculation
    [info(:,end+1) h(:,end+1) hGiven(:,end+1)] = calcInfo(data);
    
%     %Linear extrapolation
%     %Mutual information
%     fitParams = regress(mean(info,1)', [ones(length(dataFract)+1,1), 1./[round(dataFract*numOfReps)'; numOfReps]]);
%     infoExtrapolate = fitParams(1);
%     %Entropy
%     fitParams = regress(mean(h,1)', [ones(length(dataFract)+1,1), 1./[round(dataFract*numOfReps)'; numOfReps]]);
%     hExtrapolate = fitParams(1);
%     %Conditional entropy
%     fitParams = regress(mean(hGiven,1)', [ones(length(dataFract)+1,1), 1./[round(dataFract*numOfReps)'; numOfReps]]);
%     hCondExtrapolate = fitParams(1);

    %Quadratic extrapolation
    %Mutual information
    fitCurve = fit(1./[round(dataFract*numOfReps)'; numOfReps], mean(info,1)', 'poly2');
    infoExtrapolate = fitCurve.p3;
    %Entropy
    fitCurve = fit(1./[round(dataFract*numOfReps)'; numOfReps], mean(h,1)', 'poly2');
    hExtrapolate = fitCurve.p3;
    %Conditional entropy
    fitCurve = fit(1./[round(dataFract*numOfReps)'; numOfReps], mean(hGiven,1)', 'poly2');
    hCondExtrapolate = fitCurve.p3;
    
end
    
