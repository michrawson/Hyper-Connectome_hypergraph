function X_total_corr = calc_total_corr3(X,brain_region_ind,sample_ind_subset)

thresh = .00001;

samples = size(X,1);
brain_regions = size(X,2);
subjects = size(X,3);

subsamples = length(sample_ind_subset);

% X_total_corr = cell(subjects,1);
X_total_corr = zeros(   length(brain_region_ind),length(brain_region_ind),...
                        length(brain_region_ind),subjects);
parfor m = 1:subjects % parfor
    X_total_corr_m = zeros(length(brain_region_ind),length(brain_region_ind),length(brain_region_ind));
    for i_sub = 1:length(brain_region_ind) %1:brain_regions
        i = brain_region_ind(i_sub);
        for j_sub = 1:length(brain_region_ind) %1:brain_regions
            j = brain_region_ind(j_sub);
            for k_sub = 1:length(brain_region_ind) %1:brain_regions
                k = brain_region_ind(k_sub);
                for d1_sub = 1:subsamples
                    d1 = sample_ind_subset(d1_sub);
                    for d2_sub = 1:subsamples
                        d2 = sample_ind_subset(d2_sub);
                        for d3_sub = 1:subsamples
                            d3 = sample_ind_subset(d3_sub);
                        
                            agree_ind1 = abs(X(d1,i,m)-X(:,i,m))<thresh;
                            p_i = sum(double(agree_ind1))/samples;

                            agree_ind2 = abs(X(d2,j,m)-X(:,j,m))<thresh;
                            p_j = sum(double(agree_ind2))/samples;
                            
                            agree_ind3 = abs(X(d3,k,m)-X(:,k,m))<thresh;
                            p_k = sum(double(agree_ind3))/samples;
                            
                            agree_count = sum( double(agree_ind1)...
                                             .*double(agree_ind2)...
                                             .*double(agree_ind3));

                            p_ijk = agree_count/samples;
                            if p_ijk > 0
                                assert(p_i*p_j*p_k~=0);
                                X_total_corr_m(i,j,k) = X_total_corr_m(i,j,k) ...
                                    + p_ijk * log(p_ijk/(p_i*p_j*p_k));
                            end
                        end
                    end
                end
            end
        end
    end
%     X_total_corr{m,1} = X_total_corr_m;
    X_total_corr(:,:,:,m) = X_total_corr_m;
end
