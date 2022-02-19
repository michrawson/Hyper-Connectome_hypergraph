clear 
close all

for samples = [20, 40, 80]

use_subgraphs = true; % only calc for subgraph

if use_subgraphs
    load('MPRC_SZ_JimWaltz_6_8_20_timeseries_subgraph.mat', 'H', 'P_minusone');
    X_normal = H;
    X_schiz = P_minusone;

%     load('MPRC_SZ_JimWaltz_6_8_20.mat', 'H', 'P_minusone');
%     X_normal = H;
%     X_schiz = P_minusone;
%     load('MPRC_result_060820.mat','raw_idx_inlist2','raw_idx_par1','raw_idx_par2');
%     X_normal = X_normal(1:samples,raw_idx_par1,:);
%     X_schiz = X_schiz(1:samples,raw_idx_par1,:);    
else
    load('MPRC_SZ_JimWaltz_6_8_20.mat', 'H', 'P_minusone');
    X_normal = H;
    X_schiz = P_minusone;
end

X_normal = X_normal(1:samples,:,:);
X_schiz = X_schiz(1:samples,:,:);    

tic
X_normal_corr = get_corr_matrices(X_normal);
X_schiz_corr = get_corr_matrices(X_schiz);
toc

tic
X_normal_tc = get_total_corr(X_normal, X_normal_corr);
toc

tic
X_schiz_tc = get_total_corr(X_schiz,X_schiz_corr);
toc

% tic
% X_n_mi = mutual_info3(X_normal,1:brain_regions,1:samples);
% toc
% 
% tic
% X_s_mi = mutual_info3(X_schiz,1:brain_regions,1:samples);
% toc

file_name = sprintf('MPRC_SZ_JimWaltz_6_8_20_Mutual_Infor_region=raw_idx_inlist2_samples=%d_thresh=1.mat',samples);
save(file_name, 'X_normal_tc', 'X_schiz_tc');

end
