%%Generate Jensen-Shanon divergence matrix
% Input:    1. Responses in a binary raster (N_neurons x N_repeats x N_stim)
%           2. Indices of training stimuli
%           3. Indices of test stimuli
% Output:   1. Jensen-Shanon divergence matrix derived from pairwise
%           maximum entropy model
%           2. Jensen-Shanon divergence matrix derived from independent model
%           3. Hamming distance matrix

%   Copyright 2015 Elad Ganmor
function [djsMat2 djsMat1 dHamming] = genDjsMat_LargeNet(raster,trainStim,testStim)
    rand('twister',sum(100*clock)) %initialize random seed
    
    n = size(raster,1);%number of neurons
    
    trainRaster = raster(:,:,trainStim);
    testRaster = raster(:,:,testStim);

    fValsMat = genFeatureMat(n);%no need to calculate this each time
    p2RGivenS = zeros(2^n,length(trainStim));%pairwise stimulus conditioned distribution
    p1RGivenS = zeros(2^n,length(trainStim));%independent stimulus conditioned distribution
    display('Train')
    parfor s=1:length(trainStim)%go over all stimuli
        expected = getExpectedFromData(trainRaster(:,:,s));%estimate expected values (sufficient statistics for the models)
        p2RGivenS(:,s) = nesterovGradientAsscentMaxEnt3LargeNets(n,expected,fValsMat);%estimate pairwise distribution given stimulus s
        p1RGivenS(:,s) = independentModelFromExpected2(n,expected);%estimate independent distribution given stimulus s
    end
    %use bayes rule to estimate P(S|R) from P(R|S)
    p2SGivenR = transpose(bsxfun(@rdivide,p2RGivenS,sum(p2RGivenS,2)));
    p1SGivenR = transpose(bsxfun(@rdivide,p1RGivenS,sum(p1RGivenS,2)));
    clear p1RGivenS p2RGivenS %free up memory

    %%%%%generate distance matrices
    [pEmp testWords] = evalPempir(reshape(testRaster,n,[]));%get emprirical distribution and unique responses in the test set
    
    %generate Hamming distance matrix for responses in test data
    display('Hamming mat')
    dHamming = int8(genHammingMat(testWords));
  
    testWords = binWord2Dec(testWords);%convert unique respnses from binary vectors to decimal numbers
    p2SGivenR = p2SGivenR(:,testWords+1);%only keep the distribution over empirically observed responses
    p1SGivenR = p1SGivenR(:,testWords+1);

    %Calculate JS divergence matrices
    display('Ind mat')
    djsMat1 = single(DjsMat2(p1SGivenR',p1SGivenR'));
    display('Pairwise mat')
    djsMat2 = DjsMat2(p2SGivenR',p2SGivenR');
end
