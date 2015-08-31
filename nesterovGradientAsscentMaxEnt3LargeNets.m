%%% Calculate maximum entropy distribution given empirical mean of pairs
%%% and single neurons using gradient ascent (Nesterov's method).

%   Copyright 2015 Elad Ganmor
% Input:  1. n - number of neurons
%         2. empirExpected - expected values from empirical distribution
%         3. fValsMat - feature matrix generated using genFeatureMat.m
%         (optional)
%         4. initLambda - initial values for parameters (optional)
%             
% Output: 1. pModel - the pairwise maximum entropy distribution given the
%         expected values
%         2. lambda - model parameters (1 x numOfFeatures)
%         3. z = the partition function
%         4. modelExpected - Expected values of the model (same as
%         empirExpected up to convergence errors)
%         5. logarithm of the partition function



function [pModel, lambda, z, modelExpected, lnZ] = nesterovGradientAsscentMaxEnt3LargeNets(n, empirExpected, fValsMat, initLambda)
    numOfFeatures = n*(n+1)/2; z = 1; lnZ = 1;
    MIN_GRADIENT = 1e-8; %gradient value for which algorithm will stop
    RATE = 0.5;%learning rate
    MAX_ITER = 1e4;
    
    %%%%%%Allocate memory%%%%%%%%%%%%%%%%%%%%%%%%
    lambda = zeros(1, numOfFeatures,'single');
    pModel = zeros(1, 2^n,'single');
    modelExpected = zeros(1, numOfFeatures, 'single');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%Initialize values%%%%%%%%%%%%%%%%%%%%%
    zero = empirExpected==0;
    one = empirExpected>=1;
    chooseInitialLambda(nargin); %%Pick startting point (at random)
    
    if nargin<3
        stateMat = zeros(2^n, n,'single');
        featureMat = zeros(numOfFeatures, n,'single'); %%Each row is a feature
        fValsMat = zeros(numOfFeatures, 2^n,'single');
        genStateMat(); %%Matrix of all possible network states (2^n x n)
        genFeatureMat(); %%Work with features in matrix form (numOfFeatures x n)
    end
    
    %%%%Start Ascending%%%%%%%%%%%%%%%%%%%
    lambda(zero) = -100;
    gradient = evalGradient(lambda);
    gradientNorm = mean(gradient.^2);
    y = lambda;%for Neterov's method
    prevLambda = lambda;
    iter = 0;
    while (gradientNorm > MIN_GRADIENT && iter < MAX_ITER) || any(gradient(one)>1e-3)
        iter = iter + 1;
        lambda = y + RATE*evalGradient(y);%update rule
        lambda(zero) = -100;%if we have some expected values of 0 force corresponding parameters to be very small
        y = lambda + ((iter-1)/(iter+2))*(lambda - prevLambda);
        prevLambda = lambda;
        
        gradient = evalGradient(lambda); %new gradient
        gradientNorm = mean(gradient.^2);
    end
    
    
    %%%%%%%%%%%%%Sub Routines%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Choose initial values for the parameters
    function chooseInitialLambda(numOfArguments)
        if numOfArguments>3
            lambda = initLambda;
        else
            tmpExpected = empirExpected;
            tmpExpected(zero) = 1e-3;
            tmpExpected(one) = 1 - 1e-3;
            lambda(1:n) = independentModelParameters(n,tmpExpected);
        end
    end

    %Generate the feature matrix for a pairwise model: i.e. for n neurons. Each
    %column i (0 <= i <= 2^n -1) corresponds to the n bit binary representation
    %of 2^i - 1, denoted as w_i. The first n rows correspond to single bit
    %features (i.e. fValsMat(i)==1 iff w_i(i)==1). From the n+1th row and on it
    %corresponds to two bit features (i.e. if w_i(j)==1 ^ w_i(k)==1)
    function genFeatureMat()
        ind = n; %feature number for pairs
        for i=1:n
            featureMat(i,i) = 1;
            for j=i+1:n
                ind = ind + 1;
                featureMat(ind, i) = 0.5;
                featureMat(ind, j) = 0.5;
            end
        end
        fValsMat = floor(featureMat*stateMat'); %(numOfFeatures x 2^n)
        clear featureMat stateMat;
    end

    function genStateMat() %generate matrix of all network states
        for i = 0 : 2^n - 1
            ind = getOneIndices(i);
            stateMat(i+1,ind)=1;
        end
    end

    %%Get the indices of the ones in the word
    function ind = getOneIndices(x)
        word = dec2bin(x);
        ind = find(word == '1');
        ind = ind + (n - length(word));
    end

    %%calculate partition function and probabilities using matrix multiplication
    function calcDistribution(lambda) 
        energy = lambda*fValsMat;
        pModel = exp(energy); %(1 x 2^n)
        z = sum(pModel);
        lnZ = sum(energy);
        pModel = pModel/z; %(1 x 2^n)
    end

    %estimate the gradient <p_i>_data - <p_i>_model
    function g = evalGradient(p)
        calcDistribution(p);
        modelExpected = (fValsMat*pModel')'; %(1 x numOfFeatures)
        g = empirExpected - modelExpected;
    end
end
