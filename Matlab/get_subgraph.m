function X = get_subgraph(X,raw_idx_inlist2,raw_idx_par1,raw_idx_par2)

    X_g1 = upper_tri_vector_sq(X(raw_idx_inlist2,raw_idx_inlist2,:,:));
    X_g2_p1 = upper_tri_vector_sq(X(raw_idx_par1,raw_idx_par1,:,:));
    X_g2_p2 = upper_tri_vector_sq(X(raw_idx_par2,raw_idx_par2,:,:));

    X = cat(1, X_g1, X_g2_p1, X_g2_p2);
    