function X_mut_info = mutual_info2(X,m,i,j,sample_ind_subset)

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

        agree_ind1 = abs(X(d1,i,m)-X(:,i,m))<thresh;
        p_i = sum(double(agree_ind1))/samples;

        agree_ind2 = abs(X(d2,j,m)-X(:,j,m))<thresh;
        p_j = sum(double(agree_ind2))/samples;

        agree_count = sum( double(agree_ind1)...
                         .*double(agree_ind2));

        p_ij = agree_count/samples;
        X_mut_info = X_mut_info + p_ij * log(p_ij/(p_i*p_j));        
    end
end
