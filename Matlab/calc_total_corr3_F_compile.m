% gfortran -O3 -fopenmp compare_total_correlation.f03 && time ./a.out 

% mex debug:
% -v -g

mex -R2018a FFLAGS='$FFLAGS -Wall -std=gnu -fopenmp' ...
FOPTIMFLAGS='-O3' ...
FDEBUGFLAGS='' ...
calc_total_corr3_F.F90

predictor_vars = 6;
samples = 30;
subjects = 30;

[X_N, X_S] = get_simulated_data(predictor_vars,samples,subjects);

X_N_subj = X_N(:,:,1);
X_N_tc = calc_total_corr3_F(X_N_subj);
