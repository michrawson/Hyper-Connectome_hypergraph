close all
clear 

stat_sig = 10

subgraph_name_cells = {'1_of_4','2_of_4','3_of_4','4_of_4'};
% samples_v = [30];

% subgraph_name_cells = {'1_of_8','2_of_8','3_of_8','4_of_8','5_of_8','6_of_8'};

% subgraph_name_cells = { '1_of_8','2_of_8','3_of_8','4_of_8',...
%                         '5_of_8','6_of_8','7_of_8','8_of_8'};
% subgraph_name_cells = { '1_of_4'};
samples_v = [20, 30, 40];

for subgraph_name_cell = subgraph_name_cells
subgraph_name = subgraph_name_cell{1}
for thresh_cell = {'.00001'}
thresh = thresh_cell{1}
for samples = samples_v
samples

for dataset_cell = {'real'} % {'simulate'} 
    dataset = dataset_cell{1}
    
    fprintf('load dataset\n');
    tic
    
    if strcmp(dataset,'simulate')
        predictor_vars = 3;
%         samples = 30;
        subjects = 1000;
        [X_N, X_S] = get_simulated_data(predictor_vars,samples,subjects);
        
        X_N_corr = get_corr_matrices(X_N);
        X_S_corr = get_corr_matrices(X_S);
        
        X_N_corr_v = upper_tri_vector(X_N_corr);
        X_S_corr_v = upper_tri_vector(X_S_corr);     
        
        X_N_tc = get_total_corr(X_N, X_N_corr);
        X_S_tc = get_total_corr(X_S, X_S_corr);

    elseif strcmp(dataset,'real')
        [X_N,X_S,X_N_corr,X_S_corr,X_N_tc,X_S_tc] = load_real_data(samples,subgraph_name);
        
        X_N_corr_v = upper_tri_vector(X_N_corr);
        X_S_corr_v = upper_tri_vector(X_S_corr);     
        
        predictor_vars = size(X_N,2);
    else
        assert(false);    
    end

% %     X_N_corr = get_corr_matrices(X_N);
% %     X_S_corr = get_corr_matrices(X_S);
% 
%     X_N_corr_v = upper_tri_vector(X_N_corr);
%     X_S_corr_v = upper_tri_vector(X_S_corr);     
%     
%     if strcmp(dataset,'real')
% %         file_name = sprintf('data/MPRC_SZ_JimWaltz_6_8_20_Mutual_Infor_region=%s_samples=%2d_thresh=%s.mat',...
% %             subgraph_name,samples,thresh);
% %         load(file_name,'X_normal_tc','X_schiz_tc')
% %         X_N_tc = X_normal_tc;
% %         X_S_tc = X_schiz_tc;
%         X_N_tc = get_total_corr(X_N, X_N_corr);
%         X_S_tc = get_total_corr(X_S, X_S_corr);
%     elseif strcmp(dataset,'simulate')    
%         X_N_tc = get_total_corr(X_N, X_N_corr);
%         X_S_tc = get_total_corr(X_S, X_S_corr);
%     else
%         assert(false);
%     end
    toc

    subj_count = size(X_N_tc,2)+size(X_S_tc,2)
    if stat_sig==1
        rand_ind = randperm(subj_count);
    end
    
    for test_percent = [1/2]
    test_percent 
    
    for mode_cell = {'svmlinear'}%, 'svmrbf', 'svmrbfauto'} % 'tree','svmlinear','svmpolynomial'
        mode = mode_cell{1};
        l2Regularization = 0; 
        momentum = 0;
        optimizeHyperparameters = 'auto'; % Uses {'BoxConstraint','KernelScale'}
        hyperparameterOptimizationOptions = struct('Verbose',0,'ShowPlots',false,...
                                                    'UseParallel',true);
        smoothing = 0;
        use_feature_vect = true;
        is_nnet = false;

        for featureset_cell = {'adj', '3wayinfo'} % {'3wayinfo','adj','raw'}
            featureset = featureset_cell{1}

            if strcmp(featureset,'adj')
                X_N_v = X_N_corr_v;
                X_S_v = X_S_corr_v;       
            elseif strcmp(featureset,'3wayinfo')
                X_N_v = X_N_tc;
                X_S_v = X_S_tc;
            elseif strcmp(featureset,'raw')
                X_N_v = vectorize_by_subj(X_N);
                X_S_v = vectorize_by_subj(X_S);
            else
                assert(false);
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
            train_F1_v = zeros(stat_sig,1);
            test_F1_v = zeros(stat_sig,1);
            train_Rsq_v = zeros(stat_sig,1);
            test_Rsq_v = zeros(stat_sig,1);

            for stat_sig_ind = 1:stat_sig

                if stat_sig~=1
                    rand_ind = randperm(subj_count);
                end

                test_ind = rand_ind(1:floor(subj_count*test_percent));
                train_ind = rand_ind(floor(subj_count*test_percent)+1:end);

                [train_accuracy, test_accuracy, train_F1, test_F1, train_Rsq, ...
                    test_Rsq] = model_tester(mode,...
                    l2Regularization,smoothing,momentum,optimizeHyperparameters,...
                    hyperparameterOptimizationOptions,use_feature_vect,is_nnet,...
                    0,X_v,Y,test_ind,train_ind);

                train_accuracy_v(stat_sig_ind,1) = train_accuracy;
                test_accuracy_v(stat_sig_ind,1) = test_accuracy;
                train_F1_v(stat_sig_ind,1) = train_F1;
                test_F1_v(stat_sig_ind,1) = test_F1;
                train_Rsq_v(stat_sig_ind,1) = train_Rsq;
                test_Rsq_v(stat_sig_ind,1) = test_Rsq;
            end
            
            if strcmp(featureset, 'adj')
                test_accuracy_v_adj = test_accuracy_v;
            elseif strcmp(featureset, '3wayinfo')
                test_accuracy_v_3wayinfo = test_accuracy_v;                
            else
                assert(false);
            end
            
fprintf(['mean: featureset: %s, mode: %s, \n ',...
    'train_accuracy: %2.0f%s, test_accuracy: %2.0f%s \n',...
    'train_F1: %2.2f, test_F1: %2.2f, train_Rsq: %2.2f, test_Rsq: %2.2f \n\n'],...
                featureset,mode,...
                mean(train_accuracy_v),'%',mean(test_accuracy_v),'%',...
                mean(train_F1_v),mean(test_F1_v),...
                mean(train_Rsq_v),mean(test_Rsq_v));
            
fprintf(['stdev: featureset: %s, mode: %s, \n ',...
    'train_accuracy: %2.3f, test_accuracy: %2.3f \n',...
    'train_F1: %2.2f, test_F1: %2.2f, train_Rsq: %2.2f, test_Rsq: %2.2f \n\n'],...
                featureset,mode,...
                std(train_accuracy_v/100),std(test_accuracy_v/100),...
                std(train_F1_v),std(test_F1_v),...
                std(train_Rsq_v),std(test_Rsq_v));
            
% fprintf(['median: featureset: %s, mode: %s, train_accuracy: %2.0f%s, test_accuracy: %2.0f%s \n',...
%     'train_F1: %2.2f, test_F1: %2.2f, train_Rsq: %2.2f, test_Rsq: %2.2f \n'],...
%                 featureset,mode,...
%                 median(train_accuracy_v),'%',median(test_accuracy_v),'%',...
%                 median(train_F1_v),median(test_F1_v),...
%                 median(train_Rsq_v),median(test_Rsq_v));
            toc
        end
        [hypo,pval,ci,stats] = ttest2(test_accuracy_v_adj, test_accuracy_v_3wayinfo)
        fprintf('mean pval: %2.4f \n', pval);
        
    end
    end
end
end
end
end
