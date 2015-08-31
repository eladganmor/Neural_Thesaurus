%Generate the feature matrix for a pairwise model: i.e. for n neurons. Each
%column i (0 <= i <= 2^n -1) corresponds to the n bit binary representation
%of 2^i - 1, denoted as w_i. The first n rows correspond to single bit
%features (i.e. fValsMat(i)==1 iff w_i(i)==1). From the n+1th row and on it
%corresponds to two bit features (i.e. if w_i(j)==1 ^ w_i(k)==1)

%   Copyright 2015 Elad Ganmor
function fValsMat = genFeatureMat(n)
    numOfFeatures = n*(n+1)/2;
    stateMat = zeros(2^n, n);
    featureMat = zeros(numOfFeatures, n); %%Each row is a feature
    fValsMat = zeros(numOfFeatures, 2^n);
    
    genStateMat();
    genFeatureMat();
    
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
        for i = [0 : 2^n - 1]
            ind = getOneIndices(i);
            stateMat(i+1,ind)=1;
        end
    end

    function ind = getOneIndices(x) %%Get the indices of the ones in the word
        word = dec2bin(x);
        ind = find(word == '1');
        ind = ind + (n - length(word));
    end
end
