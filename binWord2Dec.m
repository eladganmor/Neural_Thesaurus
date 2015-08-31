%input: Array of binary words N_neurons x N_repeats x N_stim
%output: N_repeats x N_stim array of the corresponding decimal numbers

%   Copyright 2015 Elad Ganmor
function ind = binWord2Dec(w)
    [n numOfReps numOfStim] = size(w);
    ind = zeros(numOfReps,numOfStim);
    vec = 2.^(n-1:-1:0);
    for i=1:numOfStim
        ind(:,i) = vec*w(:,:,i);
    end
end
