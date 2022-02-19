function A_v = upper_tri_vector3D_sq(A)
% A is a nxnxn array. Vectorize upper triangle of A. Returns col.

A = squeeze(A);

n = size(A,1);

A_v = [];
for i = 1:n
    B = [];
    for j = i:n
        B = cat(1, B, squeeze(A(i,j,j:n)));
    end
    A_v = cat(1, A_v, B);
end
