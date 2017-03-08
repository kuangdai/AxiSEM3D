! simplex.f90
! created by Kuangdai on 14-May-2016 
! subroutines copied from SPECFEM. 

module simplex
    
    use iso_c_binding
    
    implicit none
    
    public :: simplex_fminsearch 
    private
    
contains
    
    pure subroutine simplex_fminsearch(tau_e, tau_s, nsls, f, nf, Qval, itercount, tolf, err) &
        bind(C, name="__simplex_MOD_simplex_fminsearch")
        
        implicit none
        
        integer, intent(in) :: nsls, nf
        double precision, intent(inout) :: tau_e(nsls) 
        double precision, intent(in) :: tau_s(nsls)
        double precision, intent(in) :: f(nf), Qval    
        integer, intent(out) :: itercount, err
        double precision, intent(out) :: tolf
        
        !Internal
        integer i,j, how
        integer, parameter :: none             = 0
        integer, parameter :: initial          = 1
        integer, parameter :: expand           = 2
        integer, parameter :: reflect          = 3
        integer, parameter :: contract_outside = 4
        integer, parameter :: contract_inside  = 5
        integer, parameter :: shrink           = 6
        
        integer :: maxiter, maxfun, n, func_evals, place(nsls + 1)
        
        double precision :: rho, chi, psi, sigma, usual_delta, zero_term_delta, tolx
        double precision :: y(nsls), v(nsls, nsls + 1), fv(nsls + 1), x(nsls), vtmp(nsls, nsls + 1)
        double precision :: xbar(nsls), xr(nsls), fxr, xe(nsls), fxe, xc(nsls), fxc, fxcc, xcc(nsls)
        
        n = nsls
        
        rho   = 1.0d0
        chi   = 2.0d0
        psi   = 0.5d0
        sigma = 0.5d0
        
        maxiter = 200 * n
        maxfun  = 200 * n
        
        tolx = 1.0e-4
        tolf = 1.0e-4
        itercount = 0
        err = 0
        
        v(:, :) = 0.0d0
        fv(:)  = 0.0d0
        
        v(:, 1) = tau_e
        x = tau_e
        
        fv(1) = feval(x, tau_s, nsls, f, nf, Qval)
        
        usual_delta = 0.05
        zero_term_delta = 0.00025
        
        do j = 1, n
            y = tau_e
            if (y(j)  /= 0.0d0) then
                y(j) = (1.0d0 + usual_delta) * y(j)
            else
                y(j) = zero_term_delta
            endif
            v(:, j + 1) = y
            x(:) = y
            fv(j + 1) = feval(x, tau_s, nsls, f, nf, Qval)
        enddo
        
        call heap_sort_local(n + 1, fv, place)
        
        do i = 1, n + 1
            vtmp(:, i) = v(:, place(i))
        enddo
        v = vtmp
        
        how = initial
        itercount = 1
        func_evals = n + 1
        
        do while (func_evals < maxfun .AND. itercount < maxiter)
            
            if (max_size_simplex(v, n) <= tolx .AND. max_value(fv, n + 1) <= tolf) then
                goto 666
            endif
            how = none
            
            ! xbar = average of the n (NOT n + 1) best points
            !     xbar = sum(v(:, 1:n), 2) / n
            xbar(:) = 0.0d0
            do i = 1, n
                do j = 1, n
                    xbar(i) = xbar(i) + v(i,j)
                enddo
                xbar(i) = xbar(i) / (n * 1.0d0)
            enddo
            xr = (1 + rho) * xbar - rho * v(:, n + 1)
            x(:) = xr
            fxr = feval(x, tau_s, nsls, f, nf, Qval)
            func_evals = func_evals + 1
            if (fxr < fv(1)) then
                ! Calculate the expansion point
                xe = (1 + rho * chi) * xbar - rho * chi * v(:, n + 1)
                x = xe
                fxe = feval(x, tau_s, nsls, f, nf, Qval)
                func_evals = func_evals + 1
                if (fxe < fxr) then
                    v(:, n + 1) = xe
                    fv(n + 1) = fxe
                    how = expand
                else
                    v(:, n + 1) = xr
                    fv(n + 1) = fxr
                    how = reflect
                endif
            else ! fv(:, 1) <= fxr
                if (fxr < fv(n)) then
                    v(:, n + 1) = xr
                    fv(n + 1) = fxr
                    how = reflect
                else ! fxr >= fv(:, n)
                    ! Perform contraction
                    if (fxr < fv(n + 1)) then
                        ! Perform an outside contraction
                        xc = (1 + psi * rho) * xbar - psi * rho * v(:, n + 1)
                        x(:) = xc
                        fxc = feval(x, tau_s, nsls, f, nf, Qval)
                        func_evals = func_evals + 1
                        if (fxc <= fxr) then
                            v(:, n + 1) = xc
                            fv(n + 1) = fxc
                            how = contract_outside
                        else
                            ! perform a shrink
                            how = shrink
                        endif
                    else
                        ! Perform an inside contraction
                        xcc = (1 - psi) * xbar + psi * v(:, n + 1)
                        x(:) = xcc
                        fxcc = feval(x, tau_s, nsls, f, nf, Qval)
                        func_evals = func_evals + 1
                        if (fxcc < fv(n + 1)) then
                            v(:, n + 1) = xcc
                            fv(n + 1) = fxcc
                            how = contract_inside
                        else
                            ! perform a shrink
                            how = shrink
                        endif
                    endif
                    if (how == shrink) then
                        do j = 2, n + 1
                            v(:, j)=v(:, 1) + sigma * (v(:, j) - v(:, 1))
                            x(:) = v(:, j)
                            fv(j) = feval(x, tau_s, nsls, f, nf, Qval)
                        enddo
                        func_evals = func_evals + n
                    endif
                endif
            endif
            
            call heap_sort_local(n + 1, fv, place)
            
            do i = 1, n + 1
                vtmp(:, i) = v(:, place(i))
            enddo
            v = vtmp
            
            itercount = itercount + 1
            
        enddo
        
        if (func_evals > maxfun) then
            err = 1
        endif
        if (itercount > maxiter) then
            err = 2
        endif
        
        666 continue
        x = v(:, 1)
        tolf = fv(1)
        tau_e = x
        
    endsubroutine simplex_fminsearch
    
    pure subroutine sub_maxwell(nf, nsls, f, tau_s, tau_e, B, A)
        
        implicit none
        
        double precision, parameter :: pi = 3.141592653589793d0
        ! Input
        integer, intent(in) :: nf, nsls
        double precision, dimension(nf), intent(in) :: f
        double precision, dimension(nsls), intent(in) :: tau_s, tau_e
        ! Output
        double precision, dimension(nf), intent(out) :: A, B
        
        integer i, j
        double precision w, demon
        
        A(:) = 1.0d0 - nsls * 1.0d0
        B(:) = 0.0d0
        do i = 1, nf
            w = 2.0d0 * pi * 10 ** f(i)
            do j = 1, nsls
                demon = 1.0d0 + w ** 2 * tau_s(j) ** 2
                A(i) = A(i) + (1.0d0 + w ** 2 * tau_e(j) * tau_s(j)) / demon
                B(i) = B(i) + w * (tau_e(j) - tau_s(j)) / demon
            enddo
        enddo
        
    endsubroutine sub_maxwell
    
    pure function feval(tau_e, tau_s, nsls, f, nf, Qval)
        
        implicit none
        ! Input
        integer, intent(in) :: nsls, nf
        double precision, dimension(nsls), intent(in)  :: tau_e, tau_s
        double precision, dimension(nf), intent(in) :: f
        double precision, intent(in) :: Qval
        double precision :: feval
        
        double precision, dimension(nf) :: A, B, tan_delta
        
        integer i
        double precision xi, iQ2, iQ
        
        call sub_maxwell(nf, nsls, f, tau_s, tau_e, B, A)
        tan_delta = B / A
        
        feval = 0.0d0
        iQ = 1.d0 / Qval
        iQ2 = (iQ) ** 2
        do i = 1, nf
            xi = sqrt((tan_delta(i) - iQ) **  2 / iQ2)
            feval = feval + xi
        enddo
        
    end function feval
    
    pure function max_value(fv, n)
        
        implicit none
        
        integer, intent(in) :: n
        double precision, intent(in) :: fv(n)
        double precision :: max_value
        integer i
        double precision m, z
        
        m = 0.0d0
        do i = 2, n
            z = abs(fv(1) - fv(i))
            if (z > m) then
                m = z
            endif
        enddo
        
        max_value = m
    end function max_value
    
    pure function max_size_simplex(v, n)
        
        implicit none
        
        integer, intent(in) :: n
        double precision, intent(in) :: v(n, n + 1)
        double precision :: max_size_simplex
        
        integer i, j
        double precision m, z
        
        m = 0.0d0
        do i = 1, n
            do j = 2, n + 1
                z = abs(v(i,j) - v(i, 1))
                if (z > m) then
                    m = z
                endif
            enddo
        enddo
        
        max_size_simplex = m
        
    end function max_size_simplex
    
    pure subroutine heap_sort_local(N, X, Y)
        
        implicit none
        
        integer, intent(in) :: N
        double precision, dimension(N), intent(inout) :: X
        integer, dimension(N), intent(out) :: Y
        
        ! local parameters
        double precision :: tmp
        integer :: itmp
        integer :: i
        
        do i = 1, n
            Y(i) = i
        enddo
        
        ! checks if anything to do
        if (N < 2) return
        
        ! builds heap
        do i = N / 2, 1, -1
            call heap_sort_siftdown_local(i, N, X, Y)
        enddo
        
        ! sorts array
        do i = N, 2, -1
            ! swaps last and first entry in this section
            tmp = X(1)
            X(1) = X(i)
            X(i) = tmp
            itmp = Y(1)
            Y(1) = Y(i)
            Y(i) = itmp
            call heap_sort_siftdown_local(1, i - 1, X, Y)
        enddo  
        
    endsubroutine heap_sort_local
    
    pure subroutine heap_sort_siftdown_local(start, bottom, X, Y)
        
        implicit none
        
        integer, intent(in) :: start, bottom
        double precision, intent(inout) :: X(:)
        integer, intent(inout) :: Y(:)
        
        ! local parameters
        integer :: i, j
        double precision :: xtmp
        integer :: ytmp
        
        i = start
        xtmp = X(i)
        ytmp = Y(i)
        
        j = 2 * i
        do while (j <= bottom)
            ! chooses larger value first in this section
            if (j < bottom) then
                if (X(j) <= X(j + 1)) j = j + 1
            endif
            
            ! checks if section already smaller than initial value
            if (X(j) < xtmp) exit
            
            X(i) = X(j)
            Y(i) = Y(j)
            i = j
            j = 2 * i
        enddo
        
        X(i) = xtmp
        Y(i) = ytmp
        
    endsubroutine heap_sort_siftdown_local
    
end module simplex

