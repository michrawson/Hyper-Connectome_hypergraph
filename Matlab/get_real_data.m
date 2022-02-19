function [X_N, X_S] = get_real_data(samples,subjects,subgraph_name)

load('data/MPRC_SZ_JimWaltz_6_8_20.mat','H','P_minusone');

X_N = H;
X_S = P_minusone;

load('data/MPRC_result_060820.mat','raw_idx_inlist2','raw_idx_par1','raw_idx_par2');

% sample_set_N = randsample(size(X_N,1),samples); % 1:samples;
% sample_set_N = 1:samples;
% sample_set_S = sample_set_N;

if strcmp(subgraph_name,'raw_idx_inlist2')
    if subjects>0 && samples>0
        X_N = X_N(1:samples,raw_idx_inlist2,1:subjects);
        X_S = X_S(1:samples,raw_idx_inlist2,1:subjects);
    elseif subjects>0
        X_N = X_N(:,raw_idx_inlist2,1:subjects);
        X_S = X_S(:,raw_idx_inlist2,1:subjects);
    elseif samples>0
        X_N = X_N(1:samples,raw_idx_inlist2,:);
        X_S = X_S(1:samples,raw_idx_inlist2,:);
    else
        X_N = X_N(:,raw_idx_inlist2,:);
        X_S = X_S(:,raw_idx_inlist2,:);
    end
    
elseif strcmp(subgraph_name,'raw_idx_par1')
    if subjects>0 && samples>0
        X_N = X_N(1:samples,raw_idx_par1,1:subjects);
        X_S = X_S(1:samples,raw_idx_par1,1:subjects);
    elseif subjects>0
        X_N = X_N(:,raw_idx_par1,1:subjects);
        X_S = X_S(:,raw_idx_par1,1:subjects);
    elseif samples>0
        X_N = X_N(1:samples,raw_idx_par1,:);
        X_S = X_S(1:samples,raw_idx_par1,:);
    else
        X_N = X_N(:,raw_idx_par1,:);
        X_S = X_S(:,raw_idx_par1,:);
    end
    
elseif strcmp(subgraph_name,'raw_idx_par2')
    if subjects>0 && samples>0
        X_N = X_N(1:samples,raw_idx_par2,1:subjects);
        X_S = X_S(1:samples,raw_idx_par2,1:subjects);
    elseif subjects>0
        X_N = X_N(:,raw_idx_par2,1:subjects);
        X_S = X_S(:,raw_idx_par2,1:subjects);
    elseif samples>0
        X_N = X_N(1:samples,raw_idx_par2,:);
        X_S = X_S(1:samples,raw_idx_par2,:);
    else
        X_N = X_N(:,raw_idx_par2,:);
        X_S = X_S(:,raw_idx_par2,:);
    end
else
    assert(false);
end
