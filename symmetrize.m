%Take input matrix A and make it symmetric 

%   Copyright 2015 Elad Ganmor
function B = symmetrize(A)
    n = length(A);
    B = zeros(n);
    for i = 1:n
        for j = i+1:n
            if not(A(i,j)==0)
                B(i,j) = A(i,j);
                B(j,i) = A(i,j);
            else
                B(i,j) = A(j,i);
                B(j,i) = A(j,i);
            end
        end
    end
end
