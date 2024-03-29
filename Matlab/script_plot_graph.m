close all
clear

subgraph_name = '1_of_8'
samples = 20

for dataset_cell = {'real'} % {'simulate'}
    dataset = dataset_cell{1}
    
    % Setup Data
    
    if strcmp(dataset,'simulate')
        predictor_vars = 3;
        %         samples = 30;
        subjects = 100;
        [X_N, X_S] = get_simulated_data(predictor_vars,samples,subjects);
        X_N_corr = get_corr_matrices(X_N);
        X_S_corr = get_corr_matrices(X_S);
        
        [~, X_N_tc] = get_total_corr(X_N, X_N_corr);
        [~, X_S_tc] = get_total_corr(X_S, X_S_corr);
    elseif strcmp(dataset,'real')
        %         samples = 80
        subjects = -1; % all
        [X_N, X_S] = get_real_data(samples,subjects,subgraph_name);
        
        X_N_corr = get_corr_matrices(X_N);
        X_S_corr = get_corr_matrices(X_S);
        
        [~, X_N_tc] = get_total_corr(X_N, X_N_corr);
        [~, X_S_tc] = get_total_corr(X_S, X_S_corr);
        
        %         [X_N,X_S,X_N_corr,X_S_corr,X_N_tc,X_S_tc] = load_real_data(samples,subgraph_name);
        predictor_vars = size(X_N,2);
    else
        assert(false);
    end
    
    X_N_tc_graph = zeros(size(X_N_corr));
    for subj_ind = 1:size(X_N_tc_graph,3)
        for i = 1:size(X_N_tc_graph,1)
            for j = 1:size(X_N_tc_graph,2)
                for k = 1:size(X_N_tc,3)
                    X_N_tc_graph(i,j,subj_ind) = X_N_tc_graph(i,j,subj_ind) ...
                        + X_N_tc(i,j,k,subj_ind);
                end
            end
        end
    end
    
    X_S_tc_graph = zeros(size(X_S_corr));
    for subj_ind = 1:size(X_S_tc_graph,3)
        for i = 1:size(X_S_tc_graph,1)
            for j = 1:size(X_S_tc_graph,2)
                for k = 1:size(X_S_tc,3)
                    X_S_tc_graph(i,j,subj_ind) = X_S_tc_graph(i,j,subj_ind) ...
                        + X_S_tc(i,j,k,subj_ind);
                end
            end
        end
    end
    
    n = size(X_N_corr,1);
    
    plot_range_N = 1 % :5
    plot_range_S = 1 % 2 % :5
    
    % Make Hypergraph file
    write_hypergraph_csv_file(X_N_tc, X_S_tc, plot_range_N, plot_range_S);
    
    % Make Plots
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set A_min and A_max
    A_min = -1;
    A_max = -1;
    for i = plot_range_N
        A = squeeze(abs(X_N_corr(:,:,i)));
        A = A - (A.*eye(n));
        if A_max < max(A,[],'all')
            A_max = max(A,[],'all');
        end
        if A_min > min(A,[],'all') || A_min == -1
            A_min = min(A,[],'all');
        end
    end
    for i = plot_range_S
        A = squeeze(abs(X_S_corr(:,:,i)));
        A = A - (A.*eye(n));
        if A_max < max(A,[],'all')
            A_max = max(A,[],'all');
        end
        if A_min > min(A,[],'all') || A_min == -1
            A_min = min(A,[],'all');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = plot_range_N
        g = graph(squeeze(abs(X_N_corr(:,:,i))),'omitselfloops','upper');
        LWidths = 5/2*g.Edges.Weight/max(g.Edges.Weight);
        figure
        plot(g,...'Layout','subspace','Dimension',31,...
            'Layout','force','Iterations',10,'WeightEffect','direct',...
            'NodeLabel',{},'LineWidth',LWidths)
        %         title('Normal Subject 2D');
        
        A = squeeze(abs(X_N_corr(:,:,i)));
        A = A - (A.*eye(n));
        
        figure
        imagesc(A, [A_min A_max])
        colorbar
    end
    for i = plot_range_S
        g = graph(squeeze(abs(X_S_corr(:,:,i))),'omitselfloops','upper');
        LWidths = 5/2*g.Edges.Weight/max(g.Edges.Weight);
        figure
        plot(g,...'Layout','subspace','Dimension',31,...
            'Layout','force','Iterations',10,'WeightEffect','direct',...
            'NodeLabel',{},'LineWidth',LWidths)
        %         title('Schitz. Subject 2D');
        
        A = squeeze(abs(X_S_corr(:,:,i)));
        A = A - (A.*eye(n));
        
        figure
        imagesc(A, [A_min A_max])
        colorbar
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set A_min and A_max
    A_min = -1;
    A_max = -1;
    for i = plot_range_N
        A = squeeze(abs(X_N_tc_graph(:,:,i)));
        A = A+A';
        if A_max < max(A,[],'all')
            A_max = max(A,[],'all');
        end
        if A_min > min(A,[],'all') || A_min == -1
            A_min = min(A,[],'all');
        end
    end
    for i = plot_range_S
        A = squeeze(abs(X_S_tc_graph(:,:,i)));
        A = A+A';
        if A_max < max(A,[],'all')
            A_max = max(A,[],'all');
        end
        if A_min > min(A,[],'all') || A_min == -1
            A_min = min(A,[],'all');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = plot_range_N
        g = graph(squeeze(abs(X_N_tc_graph(:,:,i))),'omitselfloops','upper');
        LWidths = 5/2*g.Edges.Weight/max(g.Edges.Weight);
        figure
        plot(g,...'Layout','subspace','Dimension',31,...
            'Layout','force','Iterations',10,'WeightEffect','direct',...
            'NodeLabel',{},'LineWidth',LWidths)
        %         title('Normal Subject 3D');
        
        A = squeeze(abs(X_N_tc_graph(:,:,i)));
        A = A+A';
        
        figure
        imagesc(A, [A_min A_max])
        colorbar
    end
    for i = plot_range_S
        g = graph(squeeze(abs(X_S_tc_graph(:,:,i))),'omitselfloops','upper');
        LWidths = 5/2*g.Edges.Weight/max(g.Edges.Weight);
        figure
        plot(g,...'Layout','subspace','Dimension',31,...
            'Layout','force','Iterations',10,'WeightEffect','direct',...
            'NodeLabel',{},'LineWidth',LWidths)
        %         title('Schitz. Subject 3D');
        
        A = squeeze(abs(X_S_tc_graph(:,:,i)));
        A = A+A';
        
        figure
        imagesc(A, [A_min A_max])
        colorbar
    end
end
