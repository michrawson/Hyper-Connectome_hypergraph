function [X_N, X_S] = get_real_data(samples,subjects,subgraph_name)

load('data/MPRC_SZ_JimWaltz_6_8_20.mat','H','P_minusone');

X_N = H;
X_S = P_minusone;

load('data/MPRC_result_060820.mat','raw_idx_inlist2','raw_idx_par1','raw_idx_par2');

% sample_set_N = randsample(size(X_N,1),samples); % 1:samples;
% sample_set_N = 1:samples;
% sample_set_S = sample_set_N;

q = floor(size(X_N,2)/4);
a = floor(size(X_N,2)/8);
if strcmp(subgraph_name,'1 of 2')
    subgraph = 1:(2*q);
elseif strcmp(subgraph_name,'2 of 2')
    subgraph = (2*q):size(X_N,2);
elseif strcmp(subgraph_name,'1 of 4')
    subgraph = 1:q;
elseif strcmp(subgraph_name,'2 of 4')
    subgraph = q:(2*q);
elseif strcmp(subgraph_name,'3 of 4')
    subgraph = (2*q):(3*q);
elseif strcmp(subgraph_name,'4 of 4')
    subgraph = (3*q):size(X_N,2);
elseif strcmp(subgraph_name,'1 of 8')
    subgraph = 1:a;
elseif strcmp(subgraph_name,'2 of 8')
    subgraph = a:2*a;
elseif strcmp(subgraph_name,'3 of 8')
    subgraph = 2*a:3*a;
elseif strcmp(subgraph_name,'4 of 8')
    subgraph = 3*a:4*a;
elseif strcmp(subgraph_name,'5 of 8')
    subgraph = 4*a:5*a;
elseif strcmp(subgraph_name,'6 of 8')
    subgraph = 5*a:6*a;
elseif strcmp(subgraph_name,'7 of 8')
    subgraph = 6*a:7*a;
elseif strcmp(subgraph_name,'8 of 8')
    subgraph = 7*a:size(X_N,2);
end

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
    if subjects>0 && samples>0
        X_N = X_N(1:samples,subgraph,1:subjects);
        X_S = X_S(1:samples,subgraph,1:subjects);
    elseif subjects>0
        X_N = X_N(:,subgraph,1:subjects);
        X_S = X_S(:,subgraph,1:subjects);
    elseif samples>0
        X_N = X_N(1:samples,subgraph,:);
        X_S = X_S(1:samples,subgraph,:);
    else
        X_N = X_N(:,subgraph,:);
        X_S = X_S(:,subgraph,:);
    end
end
