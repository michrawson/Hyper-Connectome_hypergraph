clear

subgraph_name_cells = { '1_of_8','2_of_8','3_of_8','4_of_8',...
                        '5_of_8','6_of_8','7_of_8','8_of_8',...
                        '1_of_4','2_of_4','3_of_4','4_of_4'};
samples_v = [20, 30, 40];
subjects = -1;

load('data/MPRC_SZ_JimWaltz_6_8_20.mat','H','P_minusone');

for subgraph_name_cell = subgraph_name_cells
    subgraph_name = subgraph_name_cell{1}
    for samples = samples_v
        samples

        X_N = H;
        X_S = P_minusone;

        q = floor(size(X_N,2)/4);
        a = floor(size(X_N,2)/8);
        if strcmp(subgraph_name,'1_of_2')
            subgraph = 1:(2*q);
        elseif strcmp(subgraph_name,'2_of_2')
            subgraph = (2*q):size(X_N,2);
        elseif strcmp(subgraph_name,'1_of_4')
            subgraph = 1:q;
        elseif strcmp(subgraph_name,'2_of_4')
            subgraph = q:(2*q);
        elseif strcmp(subgraph_name,'3_of_4')
            subgraph = (2*q):(3*q);
        elseif strcmp(subgraph_name,'4_of_4')
            subgraph = (3*q):size(X_N,2);
        elseif strcmp(subgraph_name,'1_of_8')
            subgraph = 1:a;
        elseif strcmp(subgraph_name,'2_of_8')
            subgraph = a:2*a;
        elseif strcmp(subgraph_name,'3_of_8')
            subgraph = 2*a:3*a;
        elseif strcmp(subgraph_name,'4_of_8')
            subgraph = 3*a:4*a;
        elseif strcmp(subgraph_name,'5_of_8')
            subgraph = 4*a:5*a;
        elseif strcmp(subgraph_name,'6_of_8')
            subgraph = 5*a:6*a;
        elseif strcmp(subgraph_name,'7_of_8')
            subgraph = 6*a:7*a;
        elseif strcmp(subgraph_name,'8_of_8')
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

        X_N_corr = get_corr_matrices(X_N);
        X_S_corr = get_corr_matrices(X_S);
    
        X_N_tc = get_total_corr(X_N, X_N_corr);
        X_S_tc = get_total_corr(X_S, X_S_corr);
        
        file_name = sprintf('data/data-samples=%d_subgraph_name=%s.mat',samples,subgraph_name);
        save(file_name,'X_N','X_S','X_N_corr','X_S_corr','X_N_tc','X_S_tc');

    end
end
