#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      implicit none
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

#if MX_HAS_INTERLEAVED_COMPLEX
      mwPointer mxGetDoubles
      mwPointer mxGetInt8s
#else
      mwPointer mxGetPr
#endif

      mwPointer x_ptr, y_ptr
      mwSize y_shape(4)
      mwSize copy_size
      mwPointer predictor_vars_p, samples_p, subjects_p
      integer*4 predictor_vars(1), samples(1), subjects(1)
      integer*4 mxClassIDFromClassName
      mwPointer mxCreateNumericArray

    interface
        subroutine calc_caller(x_ptr,y_ptr,samples,predictor_vars,subjects)
            mwPointer, intent(inout) :: x_ptr,y_ptr
            integer*4, intent(in) :: samples,predictor_vars,subjects
        end subroutine calc_caller
    end interface
!-----------------------------------------------------------------------
!     Check for proper number of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:compare_total_correlation:nInput',&
                                '4 inputs required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:compare_total_correlation:nOutput',&
                                'Too many output arguments.')
      endif

#if MX_HAS_INTERLEAVED_COMPLEX
    samples_p = mxGetInt8s(prhs(2))
    predictor_vars_p = mxGetInt8s(prhs(3))
    subjects_p = mxGetInt8s(prhs(4))
#else
    samples_p = mxGetPr(prhs(2))
    predictor_vars_p = mxGetPr(prhs(3))
    subjects_p = mxGetPr(prhs(4))
#endif

    copy_size = 1
    call mxCopyPtrToInteger4(samples_p,         samples,        copy_size)
    call mxCopyPtrToInteger4(predictor_vars_p,  predictor_vars, copy_size)
    call mxCopyPtrToInteger4(subjects_p,        subjects,       copy_size)

#if MX_HAS_INTERLEAVED_COMPLEX
      x_ptr = mxGetDoubles(prhs(1))
#else
      x_ptr = mxGetPr(prhs(1))
#endif

    y_shape = [predictor_vars(1), predictor_vars(1), predictor_vars(1), subjects(1)]
      plhs(1) = mxCreateNumericArray(4, y_shape, &
                                        mxClassIDFromClassName('double'), 0)

#if MX_HAS_INTERLEAVED_COMPLEX
      y_ptr = mxGetDoubles(plhs(1))
#else
      y_ptr = mxGetPr(plhs(1))
#endif

    call calc_caller(x_ptr,y_ptr,samples(1),predictor_vars(1),subjects(1))

      return
      end

!-----------------------------------------------------------------------

subroutine calc_caller(x_ptr,y_ptr,samples,predictor_vars,subjects)
    implicit none
    mwPointer, intent(inout) :: x_ptr,y_ptr
    integer*4, intent(in) :: samples,predictor_vars,subjects
      mwSize X_size, X_total_corr_len

      REAL*8, dimension(samples,predictor_vars,subjects) :: X_input
      REAL*8, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_total_corr

    interface
        function calc_total_corr3(predictor_vars,samples,subjects,X) &
                                                result(X_total_corr)
            integer*4, intent(in) :: predictor_vars,samples,subjects
            REAL*8, dimension(predictor_vars,samples,subjects), intent(in) :: X
            REAL*8, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_total_corr
        end function calc_total_corr3
    end interface

    X_size = samples*predictor_vars*subjects

      call mxCopyPtrToReal8(x_ptr,X_input,X_size)

    X_total_corr_len = predictor_vars*predictor_vars*predictor_vars*subjects

    X_total_corr = calc_total_corr3(predictor_vars,samples,subjects,X_input)

      call mxCopyReal8ToPtr(X_total_corr,y_ptr,X_total_corr_len)     

end subroutine calc_caller

function calc_total_corr3(predictor_vars,samples,subjects,X) &
                                        result(X_total_corr)
    implicit none
    integer*4, intent(in) :: predictor_vars,samples,subjects
    REAL*8, dimension(samples,predictor_vars,subjects), intent(in) :: X
    REAL*8, dimension(predictor_vars,predictor_vars,predictor_vars,subjects) :: X_total_corr
    integer :: m

    X_total_corr = 0

call omp_set_num_threads(16)

!$omp parallel
!$omp do
    do m = 1,subjects
        call calc_total_corr3_subj(predictor_vars,samples,subjects,X,X_total_corr,m)
    end do
!$omp end do
!$omp end parallel

end function calc_total_corr3

subroutine calc_total_corr3_subj(predictor_vars,samples,subjects,X,X_total_corr,m)
    implicit none
    integer*4, intent(in) :: predictor_vars,samples,subjects
    integer, intent(in) :: m
    REAL*8, dimension(samples,predictor_vars,subjects), intent(in) :: X
    REAL*8, dimension(predictor_vars,predictor_vars,predictor_vars,subjects), intent(inout) :: X_total_corr
    integer :: i,j,k,d1,d2,d3
    REAL*8, dimension(samples) :: p_i_ind, p_j_ind, p_k_ind
    REAL*8 :: p_i,p_j,p_k,p_ijk
    REAL*8 :: thresh = .01

    interface
        function mean(X) result(m)
            REAL*8, dimension(:), intent(in) :: X
            REAL*8 :: m
        end function mean
    end interface

    do i = 1,predictor_vars
        do j = i,predictor_vars
            do k = j,predictor_vars
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
                                X_total_corr(i,j,k,m) = &
                X_total_corr(i,j,k,m) + p_ijk * log(p_ijk/(p_i*p_j*p_k))
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine calc_total_corr3_subj

function get_simulated_data_indep(predictor_vars,samples,subjects) result(X)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    REAL*8, dimension(samples,predictor_vars,subjects) :: X
    
    call random_number(X)
    X = nint(X) ! rounding
    X = X*2-1
end function get_simulated_data_indep

function get_simulated_data_dep(predictor_vars,samples,subjects) result(X)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    integer :: i, j
    REAL*8, dimension(samples,predictor_vars,subjects) :: X
    
    call random_number(X) ! random (0,1)
    X = nint(X) ! rounding
    X = X*2-1
    
    do i = 1, subjects
        do j = 1, samples
            X(j,1:3,i) = [X(j,1,i)*X(j,2,i), X(j,2,i)*X(j,3,i), X(j,1,i)*X(j,3,i)]
        end do
    end do

end function get_simulated_data_dep

function mean(X) result(m)
    REAL*8, dimension(:), intent(in) :: X
    REAL*8 :: m
    
    m = sum(X)/size(X)

end function mean

function cor(X,Y) result(correlation)
    REAL*8, dimension(:), intent(in) :: X,Y
    REAL*8 :: correlation

    interface
        function mean(X) result(m)
            REAL*8, dimension(:), intent(in) :: X
            REAL*8 :: m
        end function mean
    end interface

    correlation = sum( (X - mean(X)) * (Y-mean(Y)) ) &
                    /((sum((X-mean(X))**2.0) &
                      *sum((Y-mean(Y))**2.0))**0.5)

end function cor

function get_corr_matrices(predictor_vars,samples,subjects, X) result(X_corr)
    implicit none
    integer, intent(in) :: predictor_vars, samples, subjects
    REAL*8, dimension(samples,predictor_vars,subjects), intent(in) :: X
    REAL*8, dimension(predictor_vars,predictor_vars,subjects) :: X_corr
    integer :: i, j, k
    
    interface
        function cor(X,Y) result(correlation)
            REAL*8, dimension(:), intent(in) :: X,Y
            REAL*8 :: correlation
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


