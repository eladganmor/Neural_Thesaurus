%Calculate Mutual inforamtion between stimulus and clustered response using
%two controls: Hamming distance clustering, and spike count clustering
%Input: 1. Responses in a binary raster (N_neurons x N_repeats x N_stim)
%       2. Which stimuli to use (the test stimuli)
%       3. Mapping between the decimal value of the word to its location in
%       the list (generated using mapWords.m)
%       4. Desired number of clusters (could be array)
%Output:    1. Clustered information using Hamming distance
%           2. Clustered information using spike count similarity
%           3. The number of clusters - should be close to but may differ
%           from the desired number of clusters
%           4-5. Same as 1-2 but using extrapolation to infinite data

%   Copyright 2015 Elad Ganmor

function [infoHamming infoCount numOfClustersHam infoHammingExt infoCountExt] = clusterHierControl(raster,test,map,numOfClusters)
    rand('twister',sum(100*clock))
    
    testRaster = raster(:,:,test);%Cluster only the test data
    n = size(testRaster,1);%number of neurons
    
    %%%%%%calculate Hamming distance matrix
    [~, testWords] = evalPempir(reshape(testRaster,n,[]));%get the unique responses in the data set
    dHamming = genHammingMat(testWords);%calculate the Hamming distance between all pairs of unique responses
    
    if nargin<6
        numOfClusters = [1:5:100 110:10:1000 1100:100:size(djsMat,1)];
    end
    
    refData = binWord2Dec(testRaster) + 1;%convert binary response vectors to decimal numbers
    
    %Initialize memory
    infoHamming = zeros(1,length(numOfClusters));
    infoHammingExt = zeros(1,length(numOfClusters));

    numOfClustersHam = zeros(1,length(numOfClusters));

    lHamming = linkage(squareform(dHamming),'average');%prepare linkage for clustering
    parfor k=1:length(numOfClusters)
        %%Cluster by Hamming distance
        cHam = cluster(lHamming,'maxclust', numOfClusters(k));
        numOfClustersHam(k) = length(unique(cHam));
        
        %%Calculate amount of info
        dataHam = zeros(size(refData));
        for i = 1:size(refData,1)%replace responses with their cluster assignment
            dataHam(i,:) = cHam(map(refData(i,:)));
        end
        infoHamming(k) = calcInfo(dataHam);%calculate mutual information between stimulus and response
        infoHammingExt(k) = calcInfoExtrapolate(dataHam,5);%extrapolate information to infinite data
    end
    
    %%%%%Check how much info is in spike count
    c = sum(testWords);%calculate spike count
    data = zeros(size(refData));
    for i = 1:size(refData,1)%replace each reponse with it's spike count
        data(i,:) = c(map(refData(i,:)));
    end
    infoCount = calcInfo(data);%mutual information between stimulus and spike count
    infoCountExt = calcInfoExtrapolate(data,5);%extrapolate information to infinite data
end
