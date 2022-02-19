
function X_tc_v = get_total_corr(X, X_corr, predictor_vars)
    X_tc = zeros(predictor_vars, predictor_vars, predictor_vars, size(X,3));
    parfor i = 1:size(X,3) %parfor
        X_N_subj = X(:,:,i);
        X_tc_subj = calc_total_corr3_F(X_N_subj);
        for j = 1:size(X_tc_subj,1)
            X_tc_subj(j,j,:) = X_corr(j,:,i);
            X_tc_subj(j,:,j) = X_corr(j,:,i);
            X_tc_subj(:,j,j) = X_corr(j,:,i);
        end
        X_tc(:,:,:,i) = X_tc_subj;
    end
    
    X_tc_v = [];
    for i = 1:size(X,3)
        X_tc_v = cat(2, X_tc_v, upper_tri_vector3D_sq(X_tc(:,:,:,i)));
    end
end
