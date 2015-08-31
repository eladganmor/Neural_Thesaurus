%%Calculate the Jensen-Shanon divergence between P(S|R) estimated from the
%%data and P(S|R_nearest) where R_nearest is the response most similar to R
%%based on the divergence matrix djs given as input.
% Input:    1. Responses in a binary raster (N_neurons x N_repeats x N_stim)
%           2. Divergence matrix between all nique responses in the input
%           raster
% Output:   1. Average Jensen-Shanon divergence between P(S|R) and P(S|R_nearest)
%           2. Average Jensen-Shanon divergence between P(S|R) and a flat
%           prior over S (i.e. P(S|R) = 1/|S|)
%           3. Expected error for the divergence in 1
%           4. Expected error for the divergence in 2
%           5. Standard Error of the Mean for estimates of P(S|R)

%   Copyright 2015 Elad Ganmor
function [d2 dPrior expectedError2 expectedErrorPrior sem] = decodeStim5_LargeNetWeightedV5_noInd(testRaster,djs)
    [n, numOfRepeats, numOfStim] = size(testRaster);
    [pEmp binTestWords] = evalPempir(reshape(testRaster,n,[]));%calulate empirical distribution of responses in the data
    testWords = binWord2Dec(binTestWords);%decimal representation of unique responses in the data
    %Create mapping between decimal representation of a binary word and the
    %position of that binary word in testWords
    map = zeros(2^n,1,'uint32');
    for i=1:length(testWords)
        map(testWords(i)+1) = i;
    end

    djs(isnan(djs)) = 1;%set NaNs to max value
    djs = djs + diag(inf*ones(1,size(djs,1)));%%just so we don't use a word as its own nearest neighbor
    
    %%generate true distributions
    fValsMat = genFeatureMat(n);%no need to calculate this each time
    pRGivenS = zeros(2^n,numOfStim);
    parfor s=1:numOfStim
        expected = getExpectedFromData(testRaster(:,:,s));%estimate expected values (sufficient statistics for the models)
        pRGivenS(:,s) = nesterovGradientAsscentMaxEnt3LargeNets(n,expected,fValsMat);%estimate pairwise distribution given stimulus s
    end
    pRGivenS = pRGivenS(testWords+1,:);
    pR = sum(pRGivenS,2);
    %generate pSGivenR using Bayes rule
    pSGivenR = transpose(bsxfun(@rdivide,pRGivenS,pR));
    pSGivenR(:,pR==0) = 0;

    prior = ones(1,numOfStim)/numOfStim;%Flat prior over stimuli
    d2 = zeros(1,length(testWords));
    dPrior = zeros(1,length(testWords));
    for w = 1:length(testWords)%go over all responses in the test data
        %%%Find the most similar word using the distance matrix
        %%pairwise
        [~, rHat] = min(djs(w,:));%rHat is the nearest neighbor to w
        pHat = pSGivenR(:,rHat);%P(S|rHat)
        d2(w) = DJS(pHat',pSGivenR(:,w)');
        %%Compare to prior distribution over stimuli (flat one)
        dPrior(w) = DJS(pSGivenR(:,w)',prior);
    end

    %%%Calculate average errors
    inds = min(djs)<0.25;%we will only use patterns with reasonably close neighbors
    expectedError2 = pEmp(inds)*d2(inds)' / sum(pEmp(inds));
    expectedErrorPrior = pEmp(inds)*dPrior(inds)' / sum(pEmp(inds));
    
    %%%Calculate error estimates for the values of P(S|R)
    numOfWords = numOfStim*numOfRepeats;
    errEst = bsxfun(@rdivide,pSGivenR.*(1-pSGivenR),round(pEmp*numOfWords));
    sem = bsxfun(@rdivide,sqrt(errEst),sqrt(round(pEmp*numOfWords)));
end
