function X_mut_info = mutual_info2_cond(X,m,i,j,k,sample_ind_subset)

thresh = .00001;

samples = size(X,1);
brain_regions = size(X,2);
subjects = size(X,3);

subsamples = length(sample_ind_subset);
                
X_mut_info = 0;
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
                             .*double(agree_ind3));

            p_ik = agree_count/samples;
            
            agree_count = sum(double(agree_ind2)...
                             .*double(agree_ind3));

            p_jk = agree_count/samples;

            agree_count = sum( double(agree_ind1)...
                             .*double(agree_ind2)...
                             .*double(agree_ind3));

            p_ijk = agree_count/samples;
            
            X_mut_info = X_mut_info + p_ijk * log(p_ijk*p_k/(p_ik*p_jk));       
        end
    end
end

