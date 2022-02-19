function A_v = upper_tri_vector(A)
% A is a nxnxt array. For each t, vectorize upper triangle of A.

n = size(A,1);

A_v = zeros(n*(n-1)/2, size(A,3));
for i = 1:size(A,3)
    B = [];
    for j = 1:n
        B = cat(2, B, A(j,j+1:n,i));
    end
    A_v(:,i) = B;
end
