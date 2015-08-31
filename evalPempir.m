%this functions gives the probabilities of the words observed in a given
%raster. 
%the words (used for later plotting)
% Input:    Binary raster of responses
% Output:   1. Empirical (histogram) probabilities of unique reponses in the data
%           2. The unique binary responses in the data (same order as in 1)
%           3. Counts of each response (instead of probability)
%           4. Division of the data set according to the number of spikes
%           in each response (used for plotting)
function [probs words counts parts] = evalPempir(raster)
    [n numOfSamples] = size(raster);
    
%   Copyright 2015 Elad Ganmor
    %for responses with 1-3 spikes (the vast majority of the data we use
    %complete 1D-3D arrays for direct access to the counts
    singles = zeros(1,n); %holds probabilities of events with exactly one spiking neuron
    doubles = zeros(n ,n); %holds probabilities of events with exactly two spiking neurons
    triples = zeros(n ,n, n); %holds probabilities of events with exactly three spiking neurons
    parts = zeros(1, 4);
    
    spikeNum = sum(raster);%number of spikes in each response
    c0 = length(find(spikeNum==0));%count of the all 0 (no spikes) response
    
    oneToThree = find(spikeNum > 0 & spikeNum < 4);%indices in which there are 1-3 spikes
    for i = oneToThree %go over all responses with 1-3 spikes
        num = spikeNum(i);
        spikingNeurons = find(raster(:,i))';
        switch num %increment the count of this reponse by 1 using direct access to the relevant array
            case 1
                singles(spikingNeurons) = singles(spikingNeurons) + 1;
            case 2
                doubles(subv2ind(size(doubles), spikingNeurons)) = doubles(subv2ind(size(doubles), spikingNeurons)) + 1;
            case 3
                triples(subv2ind(size(triples), spikingNeurons)) = triples(subv2ind(size(triples), spikingNeurons)) + 1;
        end
    end
    
    %for responses with more than 3 spikes (not many) we will use a
    %homemade hashtable with open hashing (standard hash functions don't work well for this
    %data since it is highly non-uniform)
    overThree = find(spikeNum>3);%indices of responses with > 3 spikes
    tmp = primes(length(overThree)*10);%find a prime number big enough for the array
    if ~isempty(tmp)
        arrLength = tmp(end);
    else
        arrLength = 0;
    end
    otherCount = zeros(1,arrLength);
    otherWords = false(n,arrLength);
    iter = 0;
    for i = overThree
        iter = iter + 1;
        spikingNeurons = find(raster(:,i));
        ind = insertInd(spikingNeurons);%find the position in the array using the hash function
        otherCount(ind) = otherCount(ind) + 1;
        otherWords(:, ind) = raster(:,i);
    end
    
    %merge all counts and reponses
    counts = zeros(1, length(find(singles)) + length(find(doubles)) + length(find(triples)) + length(find(otherCount)) + 1);
    words = false(n, length(find(singles)) + length(find(doubles)) + length(find(triples)) + length(find(otherCount)) + 1);
    counts(1) = c0;%zero spike count goes int first
    
    %add all single spike responses
    ind = 1;
    parts(1) = 1; %index where single spike words begin
    nonZero = find(singles);
    for i = nonZero
        ind = ind + 1;
        counts(ind) = singles(i);
        words(i, ind) = 1;
    end
    
    %add all 2 spike responses
    parts(2) = ind + 1; %index where double spike words begin
    nonZero = find(doubles)';
    for i = nonZero
        ind = ind + 1;
        counts(ind) = doubles(i)';
        words(ind2subv(size(doubles), i), ind) = 1;
    end
    
    %add all 3 spike responses
    parts(3) = ind + 1; %index where triple spike words begin
    nonZero = find(triples)';
    for i = nonZero
        ind = ind + 1;
        counts(ind) = triples(i);
        words(ind2subv(size(triples), i), ind) = 1;
    end

    %add all >3 spike responses
    parts(4) = ind + 1; %index where over three spike words begin
    nonZero = find(otherCount)';
    counts(ind+1:end) = double(otherCount(nonZero));
    words(:, ind+1:end) = otherWords(:, nonZero);
    probs = counts/numOfSamples;
    
    %%Subroutines%%%%%%%%
    
    %Hash function that takes the indices of the spiking neurons and
    %returns the entry in the hash table for that response
    function freeIndex = insertInd(spikes)
        w = false(n, 1);
        w(spikes) = 1;
        % to gain relatively uniform hashing we assume all neurons have
        % similar spike probability
        freeIndex = round(mod(spikes(1)*arrLength/(n-2) - 1, arrLength)) + 1;
        spikeInd = 1;
        while otherCount(freeIndex) && ~all(otherWords(:,freeIndex)==w) %search for an unused space
            spikeInd = spikeInd + 1;
            if spikeInd <= length(spikes)
                freeIndex = mod(freeIndex + spikes(spikeInd), arrLength) + 1;
            else
                freeIndex = mod(freeIndex*2, arrLength) + 1;
            end
        end
    end            
end
