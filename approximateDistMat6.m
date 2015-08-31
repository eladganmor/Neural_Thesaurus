%%This script tries to approximate (least squares) the distance matrix
%%D = dij using dij = \sum_k w^k*(wi^k - wj^k)^2
%%For each cluster we build a separate weight vector
% Input:    1. Responses in a binary raster (N_neurons x N_repeats x N_stim)
%           2. Divergence matrix between the responses estimated using
%           genDjsMat_LargeNet.m
%           3. Cluster assignment of the responses calculated using
%           clusterHier2.m
% Output:   1. weights for a 1st order linear approximation of the input
%           divergence matrix
%           2. Same as 1 but for a 2nd order model
%           3. Correlation between the input divergence matrix and the 1st
%           order approximation (per cluster)
%           4. Same as 3 but for a 2nd order model

%   Copyright 2015 Elad Ganmor
function [w w2 corr corr2] = approximateDistMat6(testRaster,djsMat,clustering)
    rand('twister',sum(100*clock))
    trainFract = 0.6;
    n = size(testRaster,1);

    %%%Pick clustering and calculate size fo clusters
    numOfClusters = max(clustering);
    clusterSize = zeros(1,numOfClusters);
    for i=1:numOfClusters
        clusterSize(i) = sum(clustering==i);
    end
    
    [~, binTestWords] = mapWords(reshape(testRaster,n,[]));
    chosenClusters = find(clusterSize>10 & clusterSize<100);%only use reasonably sized clusters

    w2 = zeros(length(chosenClusters),n*(n+1)/2); %the weight vector for a 2nd order model (single neurons + pairs)
    w = zeros(length(chosenClusters),n); %the weight vector for a 1st order model (single neurons + pairs)
    corr = zeros(length(chosenClusters),1); %correlation to ground truth for a 1st order model
    corr2 = zeros(length(chosenClusters),1); %correlation to ground truth for a 2nd order model
    dGround = cell(1,length(chosenClusters)); %ground truth distances
    d1 = cell(1,length(chosenClusters)); %1st order model approximate distances
    d2 = cell(1,length(chosenClusters)); %2nd order model approximate distances
    
    numOfTestSamples = 0;
    for clusterInd=1:length(chosenClusters)%go over all clusters
        numOfPoints = clusterSize(chosenClusters(clusterInd)).*(clusterSize(chosenClusters(clusterInd)) - 1)/2;%number of pairs
        A2 = zeros(numOfPoints,n*(n+1)/2);%%each row of this matrix is the squared difference between bits of a pair of words
        A = zeros(numOfPoints,n);
        d = zeros(numOfPoints,1);
        
        %get the responses assigned to the current cluster
        currInds = find(clustering==chosenClusters(clusterInd));
        words = binTestWords(:,currInds);
        ind = 0;
        for i=1:clusterSize(chosenClusters(clusterInd))
            for j=i+1:clusterSize(chosenClusters(clusterInd))
                ind = ind + 1;
                A2(ind,:) = calcDiff(words(:,i), words(:,j)); %(xi - yj).^2 and all 2nd order differences
                A(ind,:) = (words(:,i) - words(:,j)).^2; %(xi - yj).^2
                d(ind) = djsMat(currInds(i),currInds(j));%the "true" given distance
            end
        end
        badInds = isnan(d);
        A2(badInds,:) = [];
        A(badInds,:) = [];
        d(badInds) = [];
        
        %%Choose cross validation set
        num = length(d);
        p = randperm(num);
        train = p(1:round(trainFract*num));
        test = p(round(trainFract*num)+1:end);
        
        %solve fot the linear weights
        w(clusterInd,:) = A(train,:)\d(train);
        w2(clusterInd,:) = A2(train,:)\d(train);
        
        dGround{clusterInd} = d(test);%ground truth
        d1{clusterInd} = A(test,:)*w(clusterInd,:)';%1st order approx.
        d2{clusterInd} = A2(test,:)*w2(clusterInd,:)';%2nd order approx.
        numOfTestSamples = numOfTestSamples + length(test);
        
        %calculate correlation between ground truth and the approximations
        tmp = corrcoef(dGround{clusterInd},d1{clusterInd});
        corr(clusterInd) = tmp(2);
        tmp = corrcoef(dGround{clusterInd},d2{clusterInd});
        corr2(clusterInd) = tmp(2);
    end
    
    function delta = calcDiff(w1,w2) %difference between two words
        delta = zeros(1,n*(n+1)/2);
        delta(1:n) = (w1 - w2).^2;
        pos = n;
        for bit1 = 1:n
            for bit2 = bit1 + 1 : n
                pos = pos + 1;
                delta(pos) = (w1(bit1) - w2(bit2)).^2;
            end
        end
    end

end
