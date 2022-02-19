function X_N_corr = get_corr_matrices(X)

predictor_vars = size(X,2);
subjects = size(X,3);

X_N_corr = zeros(predictor_vars, predictor_vars, subjects);
for i = 1:subjects
    for j = 1:predictor_vars
        for k = 1:predictor_vars
            X_N_corr(j,k,i) = corr(X(:,j,i), X(:,k,i));
        end
    end
end
