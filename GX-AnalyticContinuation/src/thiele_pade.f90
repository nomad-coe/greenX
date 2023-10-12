module thiele_pade

    implicit none
    
    private
    public :: get_pade_approx_Wmn, compute_Wmn_pade
    
    contains
    
    ! **************************************************************************************************
    !> brief Gets the Pade approximant of the Wmn matrices
    !  This routine implements a modified version of the Thiele's reciprocal differences
    !  interpolation algorithm using a greedy strategy, ensuring that the ordering of the    
    !  included points minimizes the value of |P_n(x_{1+1}) - Wmn(x_{i+1})|
    !  The default Thiele interpolation is also included for conveniency
    !  o n_im -- number of purely imaginary frequencies
    !  o n_wac -- number of additional "near real" frequencies
    !  o total number of CD-WAC parameters (see warning)
    !  o x_ref -- array of total reference frequencies
    !  o y_im -- array of purely imaginary Wmn
    !  o y_wac -- array of additional Wmn
    !  o par -- array of the interpolation parameters
    !  o do_greedy -- whether to use the default greedy algorithm or the naive one
    ! **************************************************************************************************
    subroutine get_pade_approx_Wmn(n_im, n_wac, npar_in, x_ref, y_im, y_wac, par, do_greedy)
        
        implicit none 
        
        integer, intent(in)                                   :: n_im
        integer, intent(in)                                   :: n_wac
        integer, intent(in)                                   :: npar_in 
        complex(kind=8), dimension(n_im+n_wac), intent(inout) :: x_ref
        complex(kind=8), dimension(n_im), intent(in)          :: y_im
        complex(kind=8), dimension(n_wac), intent(in)         :: y_wac
        complex(kind=8), dimension(npar_in), intent(inout)    :: par
        logical, optional, intent(in)                         :: do_greedy
        
        ! Internal variables
        logical                                               :: my_do_greedy
        integer                                               :: n, i, idx, jdx, kdx, n_rem
        integer                                               :: i_par, npar
        integer, dimension(:), allocatable                    :: n_rem_idx
        real(kind=8)                                          :: deltap, pval
        complex(kind=8)                                       :: pval_in, x_in, y_in
        complex(kind=8), dimension(npar_in, npar_in)          :: g_func
        complex(kind=8), dimension(npar_in)                   :: xtmp, ytmp
        complex(kind=8), dimension(n_im+n_wac)                :: x, y
        
        ! Initialize variables
        xtmp = cmplx(0.0d0,0.0d0,kind=8)
        ytmp = cmplx(0.0d0,0.0d0,kind=8)
        par = cmplx(0.0d0,0.0d0,kind=8)
        g_func = cmplx(0.0d0,0.0d0,kind=8)
        
        ! Ensure the correct number of parameters in Thiele's interpolation
        npar = npar_in
        n = n_im + n_wac
        if(npar.ne.n) then
            npar=n
            write(use_unit,*) "* Warning! The number of parameters in Thiele's interpolation"&
            " must be equal to the number of reference points. Resetting."
        endif
        
        ! Pack pure imaginary and complex reference frequencies in a single array
        do i = 1, n_im
            y(i) = y_im(i)
        end do
        do i = 1, n_wac
            y(n_im+i) = y_wac(i)
        end do
        
        ! Whether to perform the refined Thiele's interpolation (default)
        my_do_greedy = .True.
        if (present(do_greedy)) then
            my_do_greedy = do_greedy
        end if
        
        ! Allocate auxiliary array
        if (.not.allocated(n_rem_idx)) then
            allocate(n_rem_idx(npar))
        end if
        n_rem_idx(:) = (/(i, i=1,npar)/)
        
        if (my_do_greedy) then 
            ! Unpack initial reference arguments, as they will be overwritten
            x = x_ref
            x_ref = cmplx(0.0d0,0.0d0,kind=8)
            
            ! Select first point as to maximize |Wmn|
            kdx = minloc(abs(x),dim=1)
            xtmp(1) = x(kdx)
            ytmp(1) = y(kdx)
            x_ref(1) = x(kdx)
            
            n_rem = npar - 1
            do i = kdx, n_rem 
                n_rem_idx(i) = n_rem_idx(i+1)
            end do
            
            ! Compute the generating function for the first time
            call thiele_gcoeff(xtmp,ytmp,g_func,1)
            par(1) = g_func(1,1)
            
            ! Add remaining points ensuring min |P_i(x_{1+1}) - Wmn(x_{i+1})|
            do idx = 2, npar
                pval = 6.02d23 ! A very large value to compare with
                do jdx = 1, n_rem
                    ! Compute next convergent P_i(x_{i+1})
                    call compute_Wmn_pade(idx-1,xtmp(1:idx-1),x(n_rem_idx(jdx)),idx-1,par,pval_in)
                    
                    ! Select the point that minimizes difference's absolute value
                    deltap = abs(pval_in - y(n_rem_idx(jdx)))
                    if (deltap .lt. pval) then
                        pval = deltap
                        x_in = x(n_rem_idx(jdx))
                        y_in = y(n_rem_idx(jdx))
                        kdx = jdx
                    end if
                end do
                
                ! Update indexes of non-visited points
                n_rem = n_rem - 1
                do i = kdx, n_rem 
                    n_rem_idx(i) = n_rem_idx(i+1)
                end do
                
                ! Add the winning point and recompute generating function
                x_ref(idx) = x_in
                xtmp(idx) = x_in
                ytmp(idx) = y_in
                call thiele_gcoeff(xtmp,ytmp,g_func,idx)
                
                ! Unpack parameters a_i = g_i(w_i)
                do i_par = 1, idx
                    par(i_par) = g_func(i_par,i_par)
                enddo
            end do
        else
            ! Set temporary array values
            xtmp(:) = x_ref(:)
            ytmp(:) = y(:)
            
            ! Interpolate
            call thiele_gcoeff(xtmp,ytmp,g_func,1)
            par(1) = g_func(1,1)
            do i_par = 2, npar
                call thiele_gcoeff(xtmp,ytmp,g_func,i_par)
                par(i_par) = g_func(i_par,i_par)
            enddo
        end if
        
        ! Deallocate arrays
        if (allocated(n_rem_idx)) then
            deallocate(n_rem_idx)
        end if
        
    end subroutine get_pade_approx_Wmn
    
    ! **************************************************************************************************
    !> brief Small subroutine to compute recurrence coefficients from Thiele's continued fraction
    !  This routine uses tabulation in order to efficienly compute the matrix elements g_func(:,:) 
    !  o n -- number of parameters
    !  o x -- array of the frequencies
    !  o y -- array of the Wmn matrix elements
    !  o g_func -- recurrence matrix used to compute the parameters a_n
    ! **************************************************************************************************
    subroutine thiele_gcoeff(x, y, g_func, n)
        integer, intent(in)                            :: n
        complex(kind=8), dimension(:), intent(in)      :: x, y
        complex(kind=8), dimension(:,:), intent(inout) :: g_func
        
        ! Internal variables
        integer :: idx
        
        ! Begin work (leveraging tabulation of the g_func(:,:))
        g_func(n,1) = y(n)
        if (n==1) return
        
        do idx = 2, n
            g_func(n, idx) = (g_func(idx-1,idx-1)-g_func(n,idx-1)) / &
            ((x(n)-x(idx-1))*(x(n)+x(idx-1))*g_func(n,idx-1))
        enddo
    end subroutine thiele_gcoeff
    
    ! **************************************************************************************************
    !> brief Gets the value of the Wmn matrices using the previously computed Pade approximant
    !  Here we only implement the Pade approximant evaluation     
    !  o n_freq -- is an input integer number corresponding to the number of points
    !               on the imagniary energy axis
    !  o omega -- is an input real energy containing the reference energy points 
    !  o x -- is an complex input number. 
    !  o npar -- is a positive integer input variable set to the number of parameters
    !            in the fitting function
    !  o par --  is a complex input array of length npar. On input par contains the 
    !            values of the fitting parameter
    !  o yexp --  is an output complex number, correponding to the value of the
    !             fitting function at x
    ! **************************************************************************************************
    subroutine  compute_Wmn_pade(n_freq,omega,x,npar,par,yexp)
        
        implicit none
        
        integer, intent(in)                             :: n_freq
        integer, intent(in)                             :: npar 
        complex(kind=8), intent(in)                     :: x
        complex(kind=8), dimension(n_freq), intent(in)  :: omega
        complex(kind=8), dimension(npar), intent(in)    :: par
        complex(kind=8), intent(inout)                  :: yexp
        
        ! Internal variables
        integer                                         :: i_par
        complex(kind=8)                                 :: gtmp
        
        ! Ensure right amount of points
        if(n_freq.ne.npar) then
            call aims_stop("Thiele's continued fraction needs the same number of reference frequencies&
            & as parameters!")
        endif
        
        ! Begin work
        gtmp = cmplx(1.0d0, 0.0d0, kind=8)
        
        do i_par = npar, 2, -1
            gtmp = 1.0d0 + par(i_par)*(x-omega(i_par-1))*(x+omega(i_par-1))/gtmp 
        enddo
        
        ! Compute out value
        yexp = par(1)/gtmp
        
    end subroutine compute_Wmn_pade
    
end module thiele_pade
