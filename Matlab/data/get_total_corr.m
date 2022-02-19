function X_tc_v = get_total_corr(X, X_corr)

    samples = size(X,1);
    predictor_vars = size(X,2);
    subjects = size(X,3);
    
    X_tc = zeros(predictor_vars, predictor_vars, predictor_vars, subjects);
    parfor i = 1:subjects %parfor
        X_subj = X(:,:,i);
        X_tc_subj = calc_total_corr3_subj_F(X_subj);
        for j = 1:size(X_tc_subj,1)
            X_tc_subj(j,j,:) = X_corr(j,:,i);
            X_tc_subj(j,:,j) = X_corr(j,:,i);
            X_tc_subj(:,j,j) = X_corr(j,:,i);
        end
        X_tc(:,:,:,i) = X_tc_subj;
    end
    
% FASTER: Use when memory req small 
    
%     X_tc = calc_total_corr3_F_adapter(X,samples,predictor_vars,subjects);
%     
%     parfor i = 1:subjects %parfor
%         X_tc_subj = squeeze(X_tc(:,:,:,i));
%         for j = 1:size(X_tc_subj,1)
%             X_tc_subj(j,j,:) = X_corr(j,:,i);
%             X_tc_subj(j,:,j) = X_corr(j,:,i);
%             X_tc_subj(:,j,j) = X_corr(j,:,i);
%         end
%         X_tc(:,:,:,i) = X_tc_subj;
%     end
    
%     max(abs(X_tc-X_tc2),[],'all')
    
    
    X_tc_v = [];
    for i = 1:subjects
        X_tc_v = cat(2, X_tc_v, upper_tri_vector3D_sq(X_tc(:,:,:,i)));
    end
    
end
