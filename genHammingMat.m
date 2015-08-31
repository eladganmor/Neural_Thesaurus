%Generate the Hamming distance matrix for a set of binary words in the
%input

%   Copyright 2015 Elad Ganmor
function mat = genHammingMat(words)
    [~, numOfWords] = size(words);
    mat = zeros(numOfWords,'single');
    parfor i = 1:numOfWords
        cellArr = num2cell(bsxfun(@minus,words,words(:,i)),1);%%each element in the cell array is one column to which we will apply the norm
        mat(i,:) = cellfun(@l1Norm,cellArr);
    end
    mat = symmetrize(mat);
end
