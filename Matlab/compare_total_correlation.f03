
function get_simulated_data_indep(predictor_vars,samples,subjects) result(X)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    double precision, dimension(samples,predictor_vars,subjects) :: X
    
    call random_number(X)
    X = nint(X) ! rounding
    X = X*2-1
end function get_simulated_data_indep

function get_simulated_data_dep(predictor_vars,samples,subjects) result(X)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    integer :: i, j
    double precision, dimension(samples,predictor_vars,subjects) :: X
    
    call random_number(X) ! random (0,1)
    X = nint(X) ! rounding
    X = X*2-1
    
    do i = 1, subjects
        do j = 1, samples
            X(j,:,i) = [X(j,1,i)*X(j,2,i), X(j,2,i)*X(j,3,i), X(j,1,i)*X(j,3,i), X(j,4,i), X(j,5,i), X(j,6,i)]
        end do
    end do

end function get_simulated_data_dep

function mean(X) result(m)
    double precision, dimension(:), intent(in) :: X
    double precision :: m
    
    m = sum(X)/size(X)

end function mean

function cor(X,Y) result(correlation)
    double precision, dimension(:), intent(in) :: X,Y
    double precision :: correlation

    interface
        function mean(X) result(m)
            double precision, dimension(:), intent(in) :: X
            double precision :: m
        end function mean
    end interface

    correlation = sum( (X - mean(X)) * (Y-mean(Y)) ) &
                    /((sum((X-mean(X))**2.0) &
                      *sum((Y-mean(Y))**2.0))**0.5)

end function cor

function get_corr_matrices(predictor_vars,samples,subjects, X) result(X_corr)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
    double precision, dimension(predictor_vars,predictor_vars,subjects) :: X_corr
    integer :: i, j, k
    
    interface
        function cor(X,Y) result(correlation)
            double precision, dimension(:), intent(in) :: X,Y
            double precision :: correlation
        end function cor
    end interface

    X_corr = 0

    do i = 1,subjects
        do j = 1,predictor_vars
            do k = 1,predictor_vars
                X_corr(j,k,i) = cor(X(:,j,i), X(:,k,i))
            end do
        end do
    end do

end function get_corr_matrices

function calc_total_corr3(predictor_vars,samples,subjects,m,X) result(X_total_corr)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
    double precision, dimension(predictor_vars,predictor_vars,predictor_vars) :: X_total_corr
    integer :: m,i,j,k,d1,d2,d3
    double precision, dimension(samples) :: p_i_ind, p_j_ind, p_k_ind
    double precision :: p_i,p_j,p_k,p_ijk
    double precision :: thresh = .00001

    interface
        function mean(X) result(m)
            double precision, dimension(:), intent(in) :: X
            double precision :: m
        end function mean
    end interface

    X_total_corr = 0

    do i = 1,predictor_vars
        do j = 1,predictor_vars
            do k = 1,predictor_vars
                do d1 = 1,samples

                    p_i_ind = 0
                    where (abs(X(d1,i,m)-X(:,i,m))<thresh)
                        p_i_ind = 1
                    end where
                    p_i = mean(p_i_ind)

                    do d2 = 1,samples

                        p_j_ind = 0
                        where (abs(X(d2,j,m)-X(:,j,m))<thresh)
                            p_j_ind = 1
                        end where
                        p_j = mean(p_j_ind)

                        do d3 = 1,samples

                            p_k_ind = 0
                            where (abs(X(d3,k,m)-X(:,k,m))<thresh)
                                p_k_ind = 1
                            end where
                            p_k = mean(p_k_ind)

                            p_ijk = mean(p_i_ind*p_j_ind*p_k_ind)
                            if (p_ijk > 0) then
                                X_total_corr(i,j,k) = X_total_corr(i,j,k) + p_ijk * log(p_ijk/(p_i*p_j*p_k))
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do

end function calc_total_corr3

function calc_total_corr3_subj(predictor_vars,samples,subjects,X) result(X_total_corr)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
    double precision, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_total_corr
    integer :: m

    interface
        function calc_total_corr3(predictor_vars,samples,subjects,m,X) result(X_total_corr)
            integer, intent(in) :: predictor_vars, samples, subjects
            double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
            double precision, dimension(predictor_vars,predictor_vars,predictor_vars) :: X_total_corr
        end function 
    end interface

    X_total_corr = 0

!$omp parallel
!$omp do
    do m = 1,subjects
        X_total_corr(:,:,:,m) = calc_total_corr3(predictor_vars,samples,subjects,m,X)
    end do
!$omp end do
!$omp end parallel

end function calc_total_corr3_subj

program compare_total_correlation
    implicit none
    
    interface
        function get_simulated_data_indep(predictor_vars,samples,subjects) result(X)
            integer, intent(in) :: predictor_vars, samples, subjects
            double precision, dimension(samples,predictor_vars,subjects) :: X
        end function

        function get_simulated_data_dep(predictor_vars,samples,subjects) result(X)
            integer, intent(in) :: predictor_vars, samples, subjects
            double precision, dimension(samples,predictor_vars,subjects) :: X
        end function      
        
        function get_corr_matrices(predictor_vars,samples,subjects, X) result(X_corr)
            integer, intent(in) :: predictor_vars, samples, subjects
            double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
            double precision, dimension(predictor_vars,predictor_vars,subjects) :: X_corr
        end function    

        function mean(X) result(m)
            double precision, dimension(:), intent(in) :: X
            double precision :: m
        end function       
        
        function calc_total_corr3_subj(predictor_vars,samples,subjects,X) result(X_total_corr)
            integer, intent(in) :: predictor_vars, samples, subjects
            double precision, dimension(samples,predictor_vars,subjects), intent(in) :: X
            double precision, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_total_corr
        end function 
    end interface

    integer, PARAMETER :: predictor_vars = 6
    integer, PARAMETER :: samples = 30
    integer, PARAMETER :: subjects = 30
    double precision, dimension(samples,predictor_vars,subjects) :: X_N, X_S
    double precision, dimension(predictor_vars,predictor_vars,subjects) :: X_N_corr, X_S_corr
    double precision, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_N_tc, X_S_tc
    double precision :: X_N_corr_pair, X_S_corr_pair, X_N_tc_triple, X_S_tc_triple

    X_N = get_simulated_data_dep(predictor_vars,samples,subjects)
    X_S = get_simulated_data_indep(predictor_vars,samples,subjects)
    
    X_N_corr = get_corr_matrices(predictor_vars,samples,subjects,X_N)
    X_S_corr = get_corr_matrices(predictor_vars,samples,subjects,X_S)
    
    X_N_corr_pair = mean(X_N_corr(1,2,:))
    X_S_corr_pair = mean(X_S_corr(1,2,:))
    print *, "X_N_corr_pair = ", X_N_corr_pair
    print *, "X_S_corr_pair = ", X_S_corr_pair
    
    X_N_tc = calc_total_corr3_subj(predictor_vars,samples,subjects,X_N)
    X_S_tc = calc_total_corr3_subj(predictor_vars,samples,subjects,X_S)
    
    X_N_tc_triple = mean(X_N_tc(1,2,3,:))
    X_S_tc_triple = mean(X_S_tc(1,2,3,:))
    print *, "X_N_tc_triple = ", X_N_tc_triple
    print *, "X_S_tc_triple = ", X_S_tc_triple
    
end program compare_total_correlation
