%%create a mapping from the unique words in raster (represented by it's 
% decimal value [1..2^n]) to it's position in the list of words
% Input:    1. Responses in a binary raster (N_neurons x N_samples)
% Output:   1. Decimal value of the unique words in the input
%           2. Binary representation of the unique words in the input
%           3. A mapping from the decimal representation of the word to its
%           position in binWords

%   Copyright 2015 Elad Ganmor
function [decWords binWords map] = mapWords(raster)
    n = size(raster,1);
    [~, binWords] = evalPempir(raster);
    decWords = binWord2Dec(binWords);%words we actually see in the test data
    map = zeros(2^n,1,'uint32');
    for i=1:length(decWords)
        map(decWords(i)+1) = i;
    end
end
