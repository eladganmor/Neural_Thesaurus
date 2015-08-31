%Estimate the independent model (i.e. product of marginals) for n neurons
%given the expected values calculated by getExpectedFromData.m

%   Copyright 2015 Elad Ganmor
function p1 = independentModelFromExpected2(n, expected)
    expected = expected(:);%to column shape
    p = repmat(expected(1:n),1,2^n);
    w = dec2bin(0 : 2^n - 1,n);%binary words
    silence = (w == '0')';
    p(silence) = 1 - p(silence);
    p1 = prod(p,1);%the probability of each word is the product of each neuron (bit) being 1 or 0
end
