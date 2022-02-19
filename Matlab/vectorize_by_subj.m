function X_v = vectorize_by_subj(X)

X_v = zeros(size(X,1)*size(X,2),size(X,3));
for i = 1:size(X,3)
    X_v(:,i) = reshape(X(:,:,i),[],1);
end
