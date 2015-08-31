%Use hierarchical clustering to cluster the responses in testRaster using
%the similarity matrix djsMat.
%Input: 1. Responses in a binary raster (N_neurons x N_repeats x N_stim)
%       2. Dis-Similarity matrix
%       3. Mapping between the decimal value of the word to its location in
%       the list (generated using mapWords.m)
%       4. Number of clusters we want (can be array)
%Output:    1. The cluster assignment
%           2. Mutual information between stimulus and clustered response
%           3. Mutual information between stimulus and (un-clustered) response
%           4. The number of clusters - should be close to but may differ
%           from the desired number of clusters
%           5-6. Same as 2-3, but using extrapolation to infinite data

%   Copyright 2015 Elad Ganmor
function [clusters info refInfo realNumOfClusters infoExt refInfoExt] = clusterHier2(testRaster,djsMat,map,numOfClusters)
    nanInds = all(isnan(djsMat) | djsMat==0); %Similarity matrix may contain NaNs
    numOfTestWords = length(djsMat);
    djsMat = djsMat(not(nanInds),not(nanInds));%reduce matrix only to not NaN values
    
    rand('twister',sum(100*clock))%initialize random seed

    if nargin<4
        numOfClusters = [1:5:100 110:10:1000 1100:100:size(djsMat,1)];
    end
    refData = binWord2Dec(testRaster) + 1;%convert binary raster to decimal
    refInfo = calcInfo(refData);%calculate raw mutual information between stimulus and response

    %Extrapolate information to infinite data
    info = zeros(1,length(numOfClusters));
    if nargout > 4
        infoExt = zeros(1,length(numOfClusters));
        refInfoExt = calcInfoExtrapolate(refData,4);
    end

    clusters = cell(1,length(numOfClusters));%Cell array for clusters
    realNumOfClusters = zeros(1,length(numOfClusters));

    l = linkage(squareform(djsMat),'average');%prepare linkage for clustering
    parfor k=1:length(numOfClusters)
        c = zeros(1,numOfTestWords);%cluster assignement
        cNotNaN = cluster(l,'maxclust', numOfClusters(k));%cluster only responses not associated with NaN in the similarity matrix
        
        %Put all NaN rows in the biggest cluster
        clusterNum = max(cNotNaN);
        clusterSize = zeros(1,clusterNum);
        for i = 1:clusterNum %calculate cluster size
            clusterSize(i) = sum(cNotNaN==i);
        end
        [~, biggestClust] = max(clusterSize);%get index of biggest cluster
        c(not(nanInds)) = cNotNaN;
        c(nanInds) = biggestClust;%assign NaNs to biggest cluster
        
        realNumOfClusters(k) = length(unique(c));%actual number of clusters
        clusters{k} = c;
        
        %replace responses with their cluster ID
        data = zeros(size(refData));
        for i = 1:size(refData,1)
            data(i,:) = c(map(refData(i,:)));
        end
        %mutual information between stimulus and clustered response
        info(k) = calcInfo(data);
        infoExt(k) = calcInfoExtrapolate(data,4);
    end
end
