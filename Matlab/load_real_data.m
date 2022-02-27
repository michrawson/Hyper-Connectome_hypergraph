function [X_N,X_S,X_N_corr,X_S_corr,X_N_tc,X_S_tc] = load_real_data(samples,subgraph_name)

file_name = sprintf('data/data-samples=%d_subgraph_name=%s.mat',samples,subgraph_name);
load(file_name,'X_N','X_S','X_N_corr','X_S_corr','X_N_tc','X_S_tc');
