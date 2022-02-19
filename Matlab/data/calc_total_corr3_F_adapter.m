function X_tc = calc_total_corr3_F_adapter(X,samples,predictor_vars,subjects)

predictor_vars = int8(predictor_vars);
samples = int8(samples);
subjects = int8(subjects);
X_tc = calc_total_corr3_F(X,samples,predictor_vars,subjects);
