%caluclate the parameters of an independent model (in exponential form)
%given expected values.
% Input:    1. Number of neurons
%           2. Expected value for each neuron (spike probability)
% Output:   1. Model parameters (one per neuron)
%           2. Partition function (normalization constant)
function [alpha, z] = independentModelParameters(n, expected)
    alpha = zeros(1, n);
    p = expected(1:n);
    
%   Copyright 2015 Elad Ganmor
    %z is one over the probability of all neurons being silent
    z=1;
    for i = 1 : n
        z=z*(1-p(i));
    end
    z = 1/z;

    %based on the formula p_i = exp(alpha_i)/z where p_i is the pattern
    %where all neurons are silent except the i'th neuron
    for i = 1:n
        prob = p(i); %neuron i is "on"
        for j = setdiff([1:n], i);
            prob = prob*(1-p(j));%all other neurons are off
        end
        alpha(i) = log(z*prob);
    end
end
