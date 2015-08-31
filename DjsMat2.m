% Calculate the Jensen Shannon divergence between several distributions in a matrix
% each --row-- of a matrix is a distribution.
% Input:    Matrices p,q, where each row sums to 1
% Output:   Matrix D s.t D_ij = JS(p_i || q_j), where JS stands for the
% Jensen Shanon divergence

%   Copyright 2015 Elad Ganmor
function djs = DjsMat2(p, q)
    numOfDistributions = size(p,1);
    p = p + eps;
    q = q + eps;
    djs = zeros(numOfDistributions);
    parfor i=1:numOfDistributions %calucalte the divergence between p_i and all distributions in q
        pq= bsxfun(@plus,p(i,:),q)/2; %pq_j = (p_i + q_j)/2
        d1 = p(i,:)*transpose(log2(bsxfun(@rdivide,p(i,:),pq))); %d1_j = KL(p_i || pq_j) (KL stands for Kullback-Leibler divergence)
        d2 = diag(q*transpose(log2(bsxfun(@rdivide,q,pq)))); %d2_j = KL(q_j || pq_j) (KL stands for Kullback-Leibler divergence)
        djs(i,:) = (d1 + d2')/2;
    end
end
