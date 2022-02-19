close all
clear 

stat_sig = 60
for subgraph_name_cell = {'raw_idx_inlist2'}
subgraph_name = subgraph_name_cell{1}
for thresh_cell = {'.00001'}
thresh = thresh_cell{1}
for samples = [20]
samples
for mode_cell = {'tree','svmlinear','svmrbf'} % 'svmpolynomial'
mode = mode_cell{1};
l2Regularization = 0; 
momentum = 0;
optimizeHyperparameters = 'auto'; % Uses {'BoxConstraint','KernelScale'}
hyperparameterOptimizationOptions = struct('Verbose',0,'ShowPlots',false,...
                                            'UseParallel',true);
smoothing = 0;
use_feature_vect = true;
is_nnet = false;
X = 0;

for dataset_cell = {'real'} % {'simulate','real'}
    dataset = dataset_cell{1};
    
    if strcmp(dataset,'simulate')
        predictor_vars = 40;
%         samples = 30;
        subjects = 30;
        [X_N, X_S] = get_simulated_data(predictor_vars,samples,subjects);
    elseif strcmp(dataset,'real')
%         samples = 80
        subjects = -1; % all
        [X_N, X_S] = get_real_data(samples,subjects,subgraph_name);
        predictor_vars = size(X_N,2);
    else
        assert(false);    
    end

%     X_N_corr = get_corr_matrices(X_N);
%     X_S_corr = get_corr_matrices(X_S);
% 
%     X_N_corr_v = upper_tri_vector(X_N_corr);
%     X_S_corr_v = upper_tri_vector(X_S_corr);     
               
    file_name = sprintf('data/MPRC_SZ_JimWaltz_6_8_20_Mutual_Infor_region=%s_samples=%2d_thresh=%s.mat',...
        subgraph_name,samples,thresh);
    load(file_name,'X_normal_tc','X_schiz_tc')
%     X_N_tc = get_total_corr(X_N, X_N_corr, predictor_vars);
%     X_S_tc = get_total_corr(X_S, X_S_corr, predictor_vars);
    X_N_tc = X_normal_tc;
    X_S_tc = X_schiz_tc;

    subj_count = size(X_N_tc,2)+size(X_S_tc,2);
    if stat_sig==1
        rand_ind = randperm(subj_count);
    end
    
    for featureset_cell = {'raw'} % {'3wayinfo','adj','raw'}
        featureset = featureset_cell{1};
        
        if strcmp(featureset,'adj')
            X_N_v = X_N_corr_v;
            X_S_v = X_S_corr_v;       
        elseif strcmp(featureset,'3wayinfo')
            X_N_v = X_N_tc;
            X_S_v = X_S_tc;
        elseif strcmp(featureset,'raw')
            
            X_N_v = vectorize_by_subj(X_N);
            X_S_v = vectorize_by_subj(X_S);
            
        end
        
        tic
        
        Y_N = ones(size(X_N_v,2), 1);
        Y_S = 2*ones(size(X_S_v,2), 1);

        X_v = cat(2, X_N_v, X_S_v);

        assert(subj_count == size(X_v,2));

        Y = cat(1, Y_N, Y_S);
        Y = categorical(Y);
        
        train_accuracy_v = zeros(stat_sig,1);
        test_accuracy_v = zeros(stat_sig,1);

        for stat_sig_ind = 1:stat_sig

            if stat_sig~=1
                rand_ind = randperm(subj_count);
            end
            
            test_ind = rand_ind(1:floor(subj_count/10));
            train_ind = rand_ind(floor(subj_count/10)+1:end);

            [train_accuracy, test_accuracy] = model_tester(mode,...
                l2Regularization,smoothing,momentum,optimizeHyperparameters,...
                hyperparameterOptimizationOptions,use_feature_vect,is_nnet,...
                X,X_v,Y,test_ind,train_ind);

            train_accuracy_v(stat_sig_ind,1) = train_accuracy;
            test_accuracy_v(stat_sig_ind,1) = test_accuracy;
        end
        fprintf('featureset: %s, mode: %s, train_accuracy: %2.0f%s, test_accuracy: %2.0f%s \n',...
            featureset,mode,mean(train_accuracy_v),'%',mean(test_accuracy_v),'%');
        toc
    end
end
end
end
end
end
