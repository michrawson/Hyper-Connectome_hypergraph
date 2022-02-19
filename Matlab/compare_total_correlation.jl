using LinearAlgebra, Statistics, Random #, Distributed

# addprocs(4)

function calc_total_corr3(X, brain_region_ind, sample_ind_subset)

    thresh = .00001
    
    samples = size(X,1)
    brain_regions = size(X,2)
    subjects = size(X,3)
    
    subsamples = length(sample_ind_subset)
    
    X_total_corr = zeros(length(brain_region_ind), length(brain_region_ind), length(brain_region_ind),subjects)

    # @sync @distributed 
    for m in 1:subjects 
        # X_total_corr_m = zeros(length(brain_region_ind),length(brain_region_ind),length(brain_region_ind))
        for i_sub in 1:length(brain_region_ind) 
            i = brain_region_ind[i_sub]
            for j_sub in 1:length(brain_region_ind) 
                j = brain_region_ind[j_sub]
                for k_sub in 1:length(brain_region_ind) 
                    k = brain_region_ind[k_sub]
                    for d1_sub in 1:subsamples
                        d1 = sample_ind_subset[d1_sub]
                        for d2_sub in 1:subsamples
                            d2 = sample_ind_subset[d2_sub]
                            for d3_sub in 1:subsamples
                                d3 = sample_ind_subset[d3_sub]
                            
                                agree_ind1 = abs.(X[d1,i,m].-X[:,i,m]).<thresh
                                p_i = sum(convert(Array{Float64},agree_ind1))/samples
    
                                agree_ind2 = abs.(X[d2,j,m].-X[:,j,m]).<thresh
                                p_j = sum(convert(Array{Float64},agree_ind2))/samples
                                
                                agree_ind3 = abs.(X[d3,k,m].-X[:,k,m]).<thresh
                                p_k = sum(convert(Array{Float64},agree_ind3))/samples
                                
                                agree_count = sum( convert(Array{Float64},agree_ind1).*convert(Array{Float64},agree_ind2).*convert(Array{Float64},agree_ind3))
    
                                p_ijk = agree_count/samples
                                if p_ijk > 0
                                    # assert(p_i*p_j*p_k != 0)
                                    X_total_corr[i,j,k,m] = X_total_corr[i,j,k,m] + p_ijk * log(p_ijk/(p_i*p_j*p_k))
                                end
                            end
                        end
                    end
                end
            end
        end
        # X_total_corr[:,:,:,m] = X_total_corr_m
    end
    
    return X_total_corr
end

function get_corr_matrices(X)

    predictor_vars = size(X,2)
    subjects = size(X,3)
    
    X_corr = zeros(predictor_vars, predictor_vars, subjects)
    for i in 1:subjects
        for j in 1:predictor_vars
            for k in 1:predictor_vars
                X_corr[j,k,i] = cor(X[:,j,i], X[:,k,i])
            end
        end
    end

    return X_corr
end

function get_simulated_data()

    predictor_vars = 3
    samples = 30
    subjects = 3#14
    
    X_N = rand(1:2,samples,predictor_vars,subjects)*2 .- 3
    for i in 1:subjects
        for j in 1:samples
            X_N[j,:,i] = [X_N[j,1,i]*X_N[j,2,i], X_N[j,2,i]*X_N[j,3,i], X_N[j,1,i]*X_N[j,3,i]]
        end
    end
    
    X_S = rand(1:2,samples,predictor_vars,subjects)*2 .- 3
    return (X_N, X_S)
end

@time begin
    
(X_N, X_S) = get_simulated_data()

X_full = cat(X_N, X_S, dims=3)

X_N_corr = get_corr_matrices(X_N)
X_S_corr = get_corr_matrices(X_S)

X_N_corr_pair = mean(X_N_corr[1,2,:])
X_S_corr_pair = mean(X_S_corr[1,2,:])
println("X_N_corr_pair = $X_N_corr_pair")
println("X_S_corr_pair = $X_S_corr_pair")

X_N_tc = calc_total_corr3(X_N, 1:size(X_N,2), 1:size(X_N,1))
X_S_tc = calc_total_corr3(X_S, 1:size(X_S,2), 1:size(X_S,1))

X_N_tc_triple = mean(X_N_tc[1,2,3,:])
X_S_tc_triple = mean(X_S_tc[1,2,3,:])
println("X_N_tc_triple = $X_N_tc_triple")
println("X_S_tc_triple = $X_S_tc_triple")

end
