function X_mut_info = mutual_info3(X,brain_region_ind,sample_ind_subset)


samples = size(X,1);
brain_regions = size(X,2);
subjects = size(X,3);

X_mut_info = cell(subjects,1);
parfor m = 1:subjects % parfor
    X_mut_info_m = zeros(length(brain_region_ind),length(brain_region_ind),length(brain_region_ind));
    for i_sub = 1:length(brain_region_ind) 
        i = brain_region_ind(i_sub);
        for j_sub = 1:length(brain_region_ind) 
            j = brain_region_ind(j_sub);
            for k_sub = 1:length(brain_region_ind) 
                k = brain_region_ind(k_sub);
                
                X_mut_info_m(i,j,k) = mutual_info2(X,m,i,j,sample_ind_subset) ...
                    - mutual_info2_cond(X,m,i,j,k,sample_ind_subset);
            end
        end
    end
    X_mut_info{m,1} = X_mut_info_m;
end
