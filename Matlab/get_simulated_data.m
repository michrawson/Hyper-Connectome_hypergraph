function [X_N, X_S] = get_simulated_data(predictor_vars,samples,subjects)

X_N = randi(2,samples,predictor_vars,subjects)*2-3;
for i = 1:subjects
    for j = 1:samples
        X_N(j,1:3,i) =   [X_N(j,1,i)*X_N(j,2,i)
                        X_N(j,2,i)*X_N(j,3,i)
                        X_N(j,1,i)*X_N(j,3,i)];
        if predictor_vars > 3
            X_N(j,4:end,i) = X_N(j,4:end,i)/1000000;
        end
    end
end

X_S = randi(2,samples,predictor_vars,subjects)*2-3;
