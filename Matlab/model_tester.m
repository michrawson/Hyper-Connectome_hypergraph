function [train_accuracy, test_accuracy, train_F1, test_F1, train_Rsq, test_Rsq] ...
    = model_tester(mode,...
    l2Regularization,smoothing,momentum,optimizeHyperparameters,...
    hyperparameterOptimizationOptions,use_feature_vect,is_nnet,X,X_v,Y,...
    test_ind,train_ind)

if use_feature_vect
    X_train = X_v(:,train_ind);
    X_test = X_v(:,test_ind);
else
    X_train = X(:,:,:,train_ind);
    X_test = X(:,:,:,test_ind);
end
Y_train = Y(train_ind, 1);
Y_test = Y(test_ind, 1);

% Smooth X_train
if smoothing > 0 && not(use_feature_vect)
    for i = 1:size(X_train,4)
        X_train(:,:,:,i) = X_train(:,:,:,i) ...
            + smoothing*graph_del2(squeeze(X_train(:,:,:,i)));
    end
    for i = 1:size(X_test,4)
        X_test(:,:,:,i) = X_test(:,:,:,i) ...
            + smoothing*graph_del2(squeeze(X_test(:,:,:,i)));
    end
end

if not(is_nnet)
    if strcmp(mode,'mnrfit')
        
        Mdl = mnrfit(X_train',Y_train);
        
        assert(false);
        
    elseif strcmp(mode,'svmrbf')
        
        Mdl = fitcsvm(X_train',Y_train, 'Standardize',true,...
            'KernelFunction','rbf',...
            'Solver','ISDA');
        
    elseif strcmp(mode,'svmrbfauto')
        
        Mdl = fitcsvm(X_train',Y_train, 'Standardize',true,...
            'OptimizeHyperparameters',optimizeHyperparameters,...
            'HyperparameterOptimizationOptions',hyperparameterOptimizationOptions,...
            'KernelFunction','rbf',...
            'Solver','ISDA');
        
    elseif strcmp(mode,'svmlinear')
        
        Mdl = fitcsvm(X_train',Y_train, 'Standardize',true,...
            ...                    'OptimizeHyperparameters',optimizeHyperparameters,...
            ...                    'HyperparameterOptimizationOptions',hyperparameterOptimizationOptions,...
            'KernelFunction','linear');
        ...                    'Solver','ISDA');
            
    elseif strcmp(mode,'svmpolynomial')
        
        Mdl = fitcsvm(X_train',Y_train, 'Standardize',true,...
            'OptimizeHyperparameters',optimizeHyperparameters,...
            'HyperparameterOptimizationOptions',hyperparameterOptimizationOptions,...
            'KernelFunction','polynomial',...
            'Solver','ISDA');
    elseif strcmp(mode,'tree')
        
        Mdl = fitctree(X_train',Y_train,...
            'OptimizeHyperparameters',optimizeHyperparameters,...
            'HyperparameterOptimizationOptions',hyperparameterOptimizationOptions);
    else
        Mdl = 0;
        assert(false);
    end
    
    Y_train_pred = predict(Mdl,X_train');
    Y_test_pred = predict(Mdl,X_test');
else
    if strcmp(mode,'alexnet')
        nnet = alexnet('Weights','none');
        layers = nnet;
        layers(1) = imageInputLayer([size(X,1) size(X,2) 1]);
        layers(23) = fullyConnectedLayer(1000);
        layers(24) = reluLayer;
        layers(25) = fullyConnectedLayer(2);
        layers(26) = softmaxLayer;
        layers(27) = classificationLayer;
        maxEpochs = 50;
        initialLearnRate = 1;
        miniBatchSize = size(X_train,4);
        
    elseif strcmp(mode,'standardimagenet')
        layers = [ ...
            imageInputLayer([size(X,1) size(X,2) 1])
            fullyConnectedLayer(1000)
            reluLayer
            fullyConnectedLayer(1000)
            reluLayer
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        maxEpochs = 80;
        initialLearnRate = 0.1;
        miniBatchSize = size(X_train,4);
        
    elseif strcmp(mode,'standardfeaturenet')
        
        X_train = X_train';
        X_test = X_test';
        
        layers = [ ...
            featureInputLayer(size(X_train,2))
            fullyConnectedLayer(1000)
            reluLayer
            fullyConnectedLayer(1000)
            reluLayer
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        maxEpochs = 80;
        initialLearnRate = 0.1;
        miniBatchSize = size(X_train,1);
    else
        layers = 0;
        maxEpochs = 0;
        initialLearnRate = 0;
        miniBatchSize = 0;
        assert(false);
    end
    
    options = trainingOptions('sgdm', ...
        'LearnRateSchedule','piecewise',...
        'LearnRateDropFactor',0.75,...
        'LearnRateDropPeriod',10,...
        'Momentum',momentum,...
        'L2Regularization',l2Regularization,...
        'MaxEpochs',maxEpochs,...
        'InitialLearnRate',initialLearnRate, ...
        'MiniBatchSize',miniBatchSize, ...
        'ExecutionEnvironment','gpu',...
        'Shuffle','every-epoch',...
        'Verbose',false);
    ...                'ValidationData',{X_test,Y_test}, ...
        ...                'ValidationFrequency',1,...
        ...                'Verbose',true, ...
        ...                'VerboseFrequency',10,...
        ...                'Plots','training-progress');
        
    nnet = trainNetwork(X_train,Y_train,layers,options);
    
    Y_train_pred = classify(nnet,X_train);
    Y_test_pred = classify(nnet,X_test);
end
%         train_accuracy = 100*length(Y_train_pred(Y_train_pred==Y_train))...
%                             /length(Y_train_pred);
%         test_accuracy = 100*length(Y_test_pred(Y_test_pred==Y_test))...
%                             /length(Y_test_pred);

[train_accuracy, train_F1, train_Rsq] = get_stats(Y_train_pred, Y_train);
[test_accuracy, test_F1, test_Rsq] = get_stats(Y_test_pred, Y_test);


