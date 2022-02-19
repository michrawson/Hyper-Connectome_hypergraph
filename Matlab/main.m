clear all
close all

stat_sig = 30;

% for nnet
l2Regularization = 0; 
momentum = 0;

% for svm
optimizeHyperparameters = 'auto'; % Uses {'BoxConstraint','KernelScale'}
hyperparameterOptimizationOptions = struct('Verbose',0,'ShowPlots',false,...
                                            'UseParallel',true);
mode_cell_v = {'svmlinear','svmpolynomial','svmrbf'};
% mode_cell_v = {'standardfeaturenet','standardimagenet'};%,'alexnet'};

load('data/MPRC_transformed_cor.mat', 'cor_NC_fisherZ_correct', 'cor_SZ_fisherZ_correct')

% Normal Subjects
X_N = reshape(cor_NC_fisherZ_correct, size(cor_NC_fisherZ_correct,1), ...
    size(cor_NC_fisherZ_correct,2), ...
    1, ...
    size(cor_NC_fisherZ_correct,3));

Y_N = ones(size(X_N,4), 1);

% Schiz. Subjects
X_S = reshape(cor_SZ_fisherZ_correct, size(cor_SZ_fisherZ_correct,1), ...
    size(cor_SZ_fisherZ_correct,2), ...
    1, ...
    size(cor_SZ_fisherZ_correct,3));

Y_S = 2*ones(size(X_S,4), 1);

X_full = cat(4, X_N, X_S);

Y = cat(1, Y_N, Y_S);
Y = categorical(Y);

subj_count = size(X_full,4);

for use_subgraphs = [true, false]
for use_eigs = [true, false]
use_subgraphs
use_eigs

tic
if use_eigs
    X_full_abs = abs(X_full);
end

if use_subgraphs
    load('data/MPRC_result_060820.mat','raw_idx_inlist2','raw_idx_par1','raw_idx_par2');
    X_v = get_subgraph(X_full_abs,raw_idx_inlist2,raw_idx_par1,raw_idx_par2);
    X = X_full_abs(raw_idx_inlist2,raw_idx_inlist2,:,:);
else
    X_v = upper_tri_vector_sq(X_full_abs);
    X = X_full_abs;
end

if use_eigs
    X_eig = zeros(size(X,1), size(X,2)+6, size(X,3), size(X,4));
    X_v_eig = zeros(size(X_v,1)+6*size(X,1), size(X,4));
    for i = 1:size(X,4)
        X_lap = squeeze(X(:,:,:,i)) - diag(squeeze(X(:,:,:,i))*ones(size(X,1),1));
        [V,D] = eigs(X_lap);

        X_eig(:,:,1,i) = cat(2, X(:,:,1,i), V);
        
        V_v = reshape(V,[],1);
        X_v_eig(:,i) = cat(1, X_v(:,i), V_v);
    end
    X_v = X_v_eig;

    X = X_eig;
end

toc

for mode_cell_ind = 1:length(mode_cell_v) % parfor

    mode = mode_cell_v{mode_cell_ind};

    if strcmp(mode,'mnrfit') || strcmp(mode,'svmrbf') || ...
            strcmp(mode,'svmlinear') || strcmp(mode,'svmpolynomial') ||...
            strcmp(mode,'standardfeaturenet')
        use_feature_vect = true;
    else
        use_feature_vect = false;
    end
    
    if strcmp(mode,'standardimagenet') || strcmp(mode,'standardfeaturenet') || ...
            strcmp(mode,'alexnet') 
        is_nnet = true;
    else
        is_nnet = false;
    end
    
    for smoothing = [0]%, .1, 1]
        if smoothing == 0 || (smoothing > 0 && not(use_feature_vect))
        
        tic
        
        train_accuracy_v = zeros(stat_sig,1);
        test_accuracy_v = zeros(stat_sig,1);
        if is_nnet
            parfor j = 1:stat_sig % parfor
                [train_accuracy, test_accuracy] = model_tester(mode,...
                    l2Regularization,smoothing,momentum,optimizeHyperparameters,...
                    hyperparameterOptimizationOptions,use_feature_vect,is_nnet,...
                    X,X_v,Y,subj_count);

                train_accuracy_v(j,1) = train_accuracy;
                test_accuracy_v(j,1) = test_accuracy;
            end
        else
            for j = 1:stat_sig % parfor
                [train_accuracy, test_accuracy] = model_tester(mode,...
                    l2Regularization,smoothing,momentum,optimizeHyperparameters,...
                    hyperparameterOptimizationOptions,use_feature_vect,is_nnet,...
                    X,X_v,Y,subj_count);

                train_accuracy_v(j,1) = train_accuracy;
                test_accuracy_v(j,1) = test_accuracy;
            end
        end
        
        fprintf('mode: %s, train_accuracy: %2.0f%s, test_accuracy: %2.0f%s \n',mode,...
                    mean(train_accuracy_v),'%',mean(test_accuracy_v),'%');

        toc
        end
    end
end
end
end
