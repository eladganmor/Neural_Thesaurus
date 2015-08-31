%Estimate the expected values of the features of a pairwise model: 
% Each column, i, corresponds to the n bit binary representation of i'th
% response in the data, denoted as w_i. The first n rows correspond to 
% single bit features (i.e. fValsMat(i)==1 iff w_i(i)==1). From the n+1th 
% row and on it corresponds to two bit features (i.e. if w_i(j)==1 ^ w_i(k)==1)
% Input: Binary raster of responses (0=silence, 1=spike) N_neurons x
% N_samples

%   Copyright 2015 Elad Ganmor
function expected = getExpectedFromData(data)
    [n, numOfSamples] = size(data);
    numOfFeatures = n*(n+1)/2;%n + nchoosek(n,2)
    
    expected = zeros(1, numOfFeatures);
    for i = 1:numOfSamples
        expected = expected + calcFeatures(data(:,i));
    end
    expected = expected/numOfSamples;
    
    function features = calcFeatures(x)
        features = zeros(1, numOfFeatures);
        
        features(1:n) = x;%single bit features
        
        %pairwise features
        indices = find(x);
        len = length(indices);
        for i1 = 1:len
            for i2 = i1 + 1:len
                m = indices(i1);
                l = indices(i2);
                pairInd = m*(n - 0.5*(m + 1)) + l;
                features(pairInd) = 1;
            end
        end
    end
end
