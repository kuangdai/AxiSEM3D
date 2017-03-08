! s20rts.f90
! created by Kuangdai on 17-May-2016 
! subroutines copied from SPECFEM. 

!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! S20RTS
!
! 3D mantle model S20RTS [Ritsema et al., 2010]
!
! Note that S20RTS uses transversely isotropic PREM as a background
! model, and that we use the PREM radial attenuation model when ATTENUATION is incorporated.
!
! reference:
!     J. Ritsema, A. Deuss, H.J. van Heijst and J.H. Woodhouse, 2010.
!     S20RTS: a degree-20 shear-velocity model for the mantle from new Rayleigh wave dispersion,
!     teleseismic traveltime and normal-mode splitting function measurements.
!     Geophys. J. Int., DOI: 10.1111/j.1365-246X.2010.04884.x
!--------------------------------------------------------------------------------------------------

module s20rts
    
    use iso_c_binding
    
    implicit none
    
    ! radii
    double precision :: RCMB_Par = 3480000.d0
    double precision :: RMOHO_Par = 6346600.d0
    double precision :: REARTH_Par = 6371000.d0
    
    ! three_d_mantle_model_constants
    integer, parameter :: NK_20 = 20
    integer, parameter :: NS_20 = 20
    
    ! model_s20rts_variables
    ! a = positive m  (radial, theta, phi) --> (k, l, m) (maybe other way around??)
    ! b = negative m  (radial, theta, phi) --> (k,l,-m)
    double precision, dimension(:, :, :), allocatable :: &
    S20RTS_V_dvs_a, S20RTS_V_dvs_b, S20RTS_V_dvp_a, S20RTS_V_dvp_b
    
    ! splines
    double precision, dimension(:), allocatable :: S20RTS_V_spknt
    double precision, dimension(:, :), allocatable :: S20RTS_V_qq0
    double precision, dimension(:, :, :), allocatable :: S20RTS_V_qq
    
    public :: initialize_s20rts, finalize_s20rts, perturb_s20rts 
    private
    
contains
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    subroutine initialize_s20rts(rcmb, rmoho, rearth, meta_data_p12, meta_data_s20) &
        bind(C, name="__s20rts_MOD_initialize_s20rts")
        
        implicit none
        
        double precision, intent(in) :: rcmb, rmoho, rearth
        double precision, intent(in) :: meta_data_p12(3549), meta_data_s20(9261)
        integer :: k, l, m, pos
        
        RCMB_Par = rcmb
        RMOHO_Par = rmoho
        REARTH_Par = rearth
        
        if (allocated(S20RTS_V_dvs_a)) return 
        
        allocate(S20RTS_V_dvs_a(0:NK_20, 0:NS_20, 0:NS_20), &
        S20RTS_V_dvs_b(0:NK_20, 0:NS_20, 0:NS_20), &
        S20RTS_V_dvp_a(0:NK_20, 0:NS_20, 0:NS_20), &
        S20RTS_V_dvp_b(0:NK_20, 0:NS_20, 0:NS_20), &
        S20RTS_V_spknt(NK_20 + 1), &
        S20RTS_V_qq0(NK_20 + 1, NK_20 + 1), &
        S20RTS_V_qq(3, NK_20 + 1, NK_20 + 1))
        
        pos = 1
        do k = 0, NK_20
            do l = 0, NS_20
                S20RTS_V_dvs_a(k, l, 0) = meta_data_s20(pos); pos = pos + 1
                do m = 1, l
                    S20RTS_V_dvs_a(k, l, m) = meta_data_s20(pos); pos = pos + 1
                    S20RTS_V_dvs_b(k, l, m) = meta_data_s20(pos); pos = pos + 1
                enddo
            enddo
        enddo
        
        ! P12 degree 12 P model from Ritsema
        pos = 1
        do k = 0, NK_20
            do l = 0, 12
                S20RTS_V_dvp_a(k, l, 0) = meta_data_p12(pos); pos = pos + 1
                do m = 1, l
                    S20RTS_V_dvp_a(k, l, m) = meta_data_p12(pos); pos = pos + 1
                    S20RTS_V_dvp_b(k, l, m) = meta_data_p12(pos); pos = pos + 1
                enddo
            enddo
            do l = 13, NS_20
                S20RTS_V_dvp_a(k, l, 0) = 0.0d0
                do m = 1, l
                    S20RTS_V_dvp_a(k, l, m) = 0.0d0
                    S20RTS_V_dvp_b(k, l, m) = 0.0d0
                enddo
            enddo
        enddo
        
        ! set up the splines used as radial basis functions by Ritsema
        call s20rts_splhsetup
        
    end subroutine initialize_s20rts
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    subroutine finalize_s20rts() &
        bind(C, name="__s20rts_MOD_finalize_s20rts")
                
        implicit none
        
        if (allocated(S20RTS_V_dvs_a)) then
            deallocate(S20RTS_V_dvs_a, S20RTS_V_dvs_b, & 
            S20RTS_V_dvp_a, S20RTS_V_dvp_b, &
            S20RTS_V_spknt, S20RTS_V_qq0, S20RTS_V_qq)
        endif
        
    end subroutine finalize_s20rts
    
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    ! !! PASS ABSOLUTE radius to me
    function perturb_s20rts(r_abs, theta, phi, r_center_abs, dvp, dvs) &
        bind(C, name="__s20rts_MOD_perturb_s20rts")
        
        implicit none
        
        double precision, intent(in) :: r_abs, theta, phi, r_center_abs
        double precision, intent(out) :: dvs, dvp
        logical :: perturb_s20rts
        
        ! local parameters
        double precision, parameter :: ZERO_ = 0.d0
        
        integer :: l, m, k
        double precision :: r_moho, r_cmb, xr, radius, r_center
        double precision :: dvs_alm, dvs_blm
        double precision :: dvp_alm, dvp_blm
        double precision :: radial_basis(0:NK_20)
        double precision :: sint, cost, x(2 * NS_20 + 1), dx(2 * NS_20 + 1)
        
        dvs = ZERO_
        dvp = ZERO_
        radius = r_abs / REARTH_Par
        r_moho = RMOHO_Par / REARTH_Par
        r_cmb = RCMB_Par / REARTH_Par
        r_center = r_center_abs / REARTH_Par
        if (r_center >= r_moho .or. r_center <= r_cmb) then 
            perturb_s20rts = .false.
            return
        endif
        ! element inside mantle but point outside
        if (radius >= r_moho * 0.999999d0) radius = r_moho * 0.999999d0 
        if (radius <= r_cmb * 1.000001d0) radius = r_cmb * 1.000001d0
        
        xr = -1.0d0 + 2.0d0 * (radius - r_cmb) / (r_moho - r_cmb)
        do k = 0, NK_20
            radial_basis(k) = s20rts_rsple(1, NK_20 + 1, S20RTS_V_spknt(1), &
            S20RTS_V_qq0(1, NK_20 + 1 - k), S20RTS_V_qq(1, 1, NK_20 + 1 - k), xr)
        enddo
        
        do l = 0, NS_20
            sint = dsin(theta)
            cost = dcos(theta)
            call lgndr(l, cost, sint, x, dx)
            
            dvs_alm = 0.0d0
            dvp_alm = 0.0d0
            do k = 0, NK_20
                dvs_alm = dvs_alm + radial_basis(k) * S20RTS_V_dvs_a(k, l, 0)
                dvp_alm = dvp_alm + radial_basis(k) * S20RTS_V_dvp_a(k, l, 0)
            enddo
            dvs = dvs + dvs_alm * x(1)
            dvp = dvp + dvp_alm * x(1)
            
            do m = 1, l
                dvs_alm = 0.0d0
                dvp_alm = 0.0d0
                dvs_blm = 0.0d0
                dvp_blm = 0.0d0
                do k = 0, NK_20
                    dvs_alm = dvs_alm + radial_basis(k) * S20RTS_V_dvs_a(k, l, m)
                    dvp_alm = dvp_alm + radial_basis(k) * S20RTS_V_dvp_a(k, l, m)
                    dvs_blm = dvs_blm + radial_basis(k) * S20RTS_V_dvs_b(k, l, m)
                    dvp_blm = dvp_blm + radial_basis(k) * S20RTS_V_dvp_b(k, l, m)
                enddo
                dvs = dvs + (dvs_alm * dcos(dble(m) * phi) + dvs_blm * dsin(dble(m) * phi)) * x(m + 1)
                dvp = dvp + (dvp_alm * dcos(dble(m) * phi) + dvp_blm * dsin(dble(m) * phi)) * x(m + 1)
            enddo
        enddo
        
        perturb_s20rts = .true.
        
    end function perturb_s20rts
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    subroutine s20rts_splhsetup
        
        implicit none
        
        ! local parameters
        integer :: i, j
        double precision :: qqwk(3, NK_20 + 1)
        
        S20RTS_V_spknt(1) = -1.00000d0
        S20RTS_V_spknt(2) = -0.78631d0
        S20RTS_V_spknt(3) = -0.59207d0
        S20RTS_V_spknt(4) = -0.41550d0
        S20RTS_V_spknt(5) = -0.25499d0
        S20RTS_V_spknt(6) = -0.10909d0
        S20RTS_V_spknt(7) = 0.02353d0
        S20RTS_V_spknt(8) = 0.14409d0
        S20RTS_V_spknt(9) = 0.25367d0
        S20RTS_V_spknt(10) = 0.35329d0
        S20RTS_V_spknt(11) = 0.44384d0
        S20RTS_V_spknt(12) = 0.52615d0
        S20RTS_V_spknt(13) = 0.60097d0
        S20RTS_V_spknt(14) = 0.66899d0
        S20RTS_V_spknt(15) = 0.73081d0
        S20RTS_V_spknt(16) = 0.78701d0
        S20RTS_V_spknt(17) = 0.83810d0
        S20RTS_V_spknt(18) = 0.88454d0
        S20RTS_V_spknt(19) = 0.92675d0
        S20RTS_V_spknt(20) = 0.96512d0
        S20RTS_V_spknt(21) = 1.00000d0
        
        do i = 1, NK_20+1
            do j = 1, NK_20+1
                if (i == j) then
                    S20RTS_V_qq0(j, i) = 1.0d0
                else
                    S20RTS_V_qq0(j, i) = 0.0d0
                endif
            enddo
        enddo
        do i = 1, NK_20 + 1
            call s20rts_rspln(1, NK_20 + 1, S20RTS_V_spknt(1), &
            S20RTS_V_qq0(1, i), S20RTS_V_qq(1, 1, i), qqwk(1, 1))
        enddo
        
    end subroutine s20rts_splhsetup
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    ! changed the obsolescent f77 features in the two routines below
    ! now still awful Fortran, but at least conforms to f90 standard
    
    double precision function s20rts_rsple(I1, I2, X, Y, Q, S)
        
        implicit none
        ! rsple returns the value of the function y(x) evaluated at point S
        ! using the cubic spline coefficients computed by rspln and saved in Q.
        ! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
        ! using the first or last interpolation polynomial. The arrays must
        ! be dimensioned at least - x(i2), y(i2), and q(3,i2).
        
        integer i1, i2
        double precision X(*), Y(*), Q(3, *), S
        
        integer i, ii
        double precision h
        
        i = 1
        II = I2 - 1
        
        !   GUARANTEE I WITHIN BOUNDS.
        I = MAX0(I, I1)
        I = MIN0(I, II)
        
        !   SEE IF X IS INCREASING OR DECREASING.
        if (X(I2) - X(I1) <  0) goto 1
        if (X(I2) - X(I1) >= 0) goto 2
        
        !   X IS DECREASING.  CHANGE I AS NECESSARY.
        1   if (S - X(I) <= 0) goto 3
        if (S - X(I) >  0) goto 4
        
        4   I = I - 1
        
        if (I - I1 <  0) goto 11
        if (I - I1 == 0) goto 6
        if (I - I1 >  0) goto 1
        
        3   if (S - X(I + 1) <  0) goto 5
        if (S - X(I + 1) >= 0) goto 6
        
        5   I = I + 1
        
        if (I - II <  0) goto 3
        if (I - II == 0) goto 6
        if (I - II >  0) goto 7
        
        !   X IS INCREASING.  CHANGE I AS NECESSARY.
        2   if (S - X(I + 1) <= 0) goto 8
        if (S - X(I + 1) >  0) goto 9
        
        9   I = I + 1
        
        if (I - II <  0) goto 2
        if (I - II == 0) goto 6
        if (I - II >  0) goto 7
        
        8   if (S - X(I) <  0) goto 10
        if (S - X(I) >= 0) goto 6
        
        10  I = I - 1
        if (I - I1 <  0) goto 11
        if (I - I1 == 0) goto 6
        if (I - I1 >  0) goto 8
        
        7   I = II
        GOTO 6
        11  I = I1
        
        !   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
        6   H = S - X(I)
        S20RTS_RSPLE = Y(I) + H * (Q(1, I) + H* (Q(2, I) + H * Q(3, I)))
        
    end function s20rts_rsple
    
    !
    !--------------------------------------------------------------------------------------------------
    !
    
    subroutine s20rts_rspln(I1, I2, X, Y, Q, F)
        
        implicit none
        ! Subroutine rspln computes cubic spline interpolation coefficients
        ! for y(x) between grid points i1 and i2 saving them in q.The
        ! interpolation is continuous with continuous first and second
        ! derivatives. It agrees exactly with y at grid points and with the
        ! three point first derivatives at both end points (i1 and i2).
        ! X must be monotonic but if two successive values of x are equal
        ! a discontinuity is assumed and separate interpolation is done on
        ! each strictly monotonic segment. The arrays must be dimensioned at
        ! least - x(i2), y(i2), q(3,i2), and f(3,i2).
        ! F is working storage for rspln.
        
        integer i1,i2
        double precision X(*), Y(*), Q(3, *), F(3, *)
        
        integer i, j, k, j1, j2
        double precision y0, a0, b0, b1, h, h2, ha, h2a, h3a, h2b
        double precision YY(3), small
        equivalence (YY(1), Y0)
        data SMALL/1.0d-08/, YY/0.0d0,0.0d0,0.0d0/
        
        J1 = I1 + 1
        Y0 = 0.0d0
        
        !   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
        if (I2 - I1  < 0) return
        if (I2 - I1 == 0) goto 17
        if (I2 - I1  > 0) goto 8
        
        8   A0 = X(J1 - 1)
        !   SEARCH FOR DISCONTINUITIES.
        DO 3 I = J1, I2
            B0 = A0
            A0 = X(I)
            if (DABS((A0 - B0) / DMAX1(A0, B0)) < SMALL) GOTO 4
            3   CONTINUE
            17  J1 = J1 - 1
            J2 = I2 - 2
            GOTO 5
            4   J1 = J1 - 1
            J2 = I - 3
            !   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
            5   if (J2 + 1 - J1 <  0) goto 9
            if (J2 + 1 - J1 == 0) goto 10
            if (J2 + 1 - J1 >  0) goto 11
            
            !   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
            10  J2 = J2 + 2
            Y0 = (Y(J2) - Y(J1)) / (X(J2) - X(J1))
            DO J=1, 3
                Q(J, J1) = YY(J)
                Q(J, J2) = YY(J)
            enddo
            GOTO 12
            
            !   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
            11  A0 = 0.
            H = X(J1 + 1) - X(J1)
            H2 = X(J1 + 2) - X(J1)
            Y0 = H * H2 * (H2 - H)
            H = H * H
            H2 = H2 * H2
            !   CALCULATE DERIVATIVE AT NEAR END.
            B0 = (Y(J1) * (H - H2) + Y(J1 + 1) * H2 - Y(J1 + 2) * H) / Y0
            B1 = B0
            
            !   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
            DO I=J1, J2
                H = X(I+1) - X(I)
                Y0 = Y(I+1) - Y(I)
                H2 = H * H
                HA = H - A0
                H2A = H - 2.0d0 * A0
                H3A = 2.0d0 * H - 3.0d0 * A0
                H2B = H2 * B0
                Q(1, I) = H2 / HA
                Q(2, I) = -HA / (H2A * H2)
                Q(3, I) = -H * H2A / H3A
                F(1, I) = (Y0 - H * B0) / (H * HA)
                F(2, I) = (H2B - Y0 * (2.0d0 * H - A0))/(H * H2 * H2A)
                F(3, I) = -(H2B - 3.0d0 * Y0 * HA) / (H * H3A)
                A0 = Q(3, I)
                B0 = F(3, I)
            enddo
            
            !   TAKE CARE OF LAST TWO ROWS.
            I = J2 + 1
            H = X(I + 1) - X(I)
            Y0 = Y(I + 1) - Y(I)
            H2 = H * H
            HA = H - A0
            H2A = H * HA
            H2B = H2 * B0 - Y0 * (2.0d0 * H - A0)
            Q(1, I) = H2 / HA
            F(1, I) = (Y0 - H * B0) / H2A
            HA = X(J2) - X(I+1)
            Y0 = -H * HA * (HA + H)
            HA = HA * HA
            
            !   CALCULATE DERIVATIVE AT FAR END.
            Y0 = (Y(I + 1) * (H2 - HA) + Y(I) * HA - Y(J2) * H2) / Y0
            Q(3, I) = (Y0 * H2A + H2B) / (H * H2 * (H - 2.0d0 * A0))
            Q(2, I) = F(1, I) - Q(1, I) * Q(3, I)
            
            !   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
            DO J=J1, J2
                K = I - 1
                Q(1, I) = F(3, K) - Q(3, K) * Q(2, I)
                Q(3, K) = F(2, K) - Q(2, K) * Q(1, I)
                Q(2, K) = F(1, K) - Q(1, K) * Q(3, K)
                I = K
            enddo
            Q(1, I) = B1
            !   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
            9   J2 = J2 + 2
            DO J = 1, 3
                Q(J, J2) = YY(J)
            enddo
            
            !   SEE IF THIS DISCONTINUITY IS THE LAST.
            12  if (J2 - I2 < 0) then
            goto 6
        else
            return
        endif
        
        !   NO.  GO BACK FOR MORE.
        6   J1 = J2 + 2
        if (J1 - I2 <= 0) goto 8
        if (J1 - I2 >  0) goto 7
        
        !   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
        7   DO J=1, 3
        Q(J, I2) = YY(J)
    enddo
    
end subroutine s20rts_rspln

!
!--------------------------------------------------------------------------------------------------
!

subroutine lgndr(l, c, s, x, dx)
    
    ! computes Legendre function x(l, m, theta)
    ! theta = colatitude, c = cos(theta), s = sin(theta), l = angular order, 
    ! sin(theta) restricted so that sin(theta) > 1.e-7
    ! x(1) contains m = 0,  x(2) contains m = 1,  x(k + 1) contains m = k
    ! m = azimuthal (longitudinal) order 0 < =  m < =  l
    ! dx = dx/dtheta
    !
    ! subroutine originally came from Physics Dept. Princeton through
    ! Peter Davis,  modified by Jeffrey Park
    
    implicit none
    
    ! argument variables
    integer l
    double precision x(2 * l + 1), dx(2 * l + 1)
    double precision c, s
    
    ! local variables
    integer i, lp1, lpsafe, lsave
    integer m, maxsin, mmm, mp1
    
    double precision sqroot2over2, c1, c2, cot
    double precision ct, d, f1, f2
    double precision f3, fac, g1, g2
    double precision g3, rfpi, sqroot3, sos
    double precision ss, stom, t, tol
    double precision v, y
    
    tol = 1.d-05
    rfpi = 0.282094791773880d0
    sqroot3 = 1.73205080756890d0
    sqroot2over2 = 0.707106781186550d0
    
    if (s >= 1.0d0 - tol) s = 1.0d0 - tol
    lsave = l
    if (l < 0) l = -1 - l
    if (l > 0) goto 1
    x(1) = rfpi
    dx(1) = 0.0d0
    l = lsave
    return
    1   if (l /=  1) goto 2
    c1 = sqroot3 * rfpi
    c2 = sqroot2over2 * c1
    x(1) = c1 * c
    x(2) = -c2 * s
    dx(1) = -c1 * s
    dx(2) = -c2 * c
    l = lsave
    return
    2   sos = s
    if (s < tol) s = tol
    cot = c / s
    ct = 2.0d0 * c
    ss = s * s
    lp1 = l + 1
    g3 = 0.0d0
    g2 = 1.0d0
    f3 = 0.0d0
    
    ! evaluate m = l value,  sans (sin(theta)) *  * l
    do i = 1, l
        g2 = g2 * (1.0d0 - 1.0d0 / (2.0d0 * i))
    enddo
    g2 = rfpi * dsqrt((2 * l + 1) * g2)
    f2 = l * cot * g2
    x(lp1) = g2
    dx(lp1) = f2
    v = 1.0d0
    y = 2.0d0 * l
    d = dsqrt(v * y)
    t = 0.0d0
    mp1 = l
    m = l - 1
    
    ! these recursions are similar to ordinary m-recursions,  but since we
    ! have taken the s *  * m factor out of the xlm's,  the recursion has the powers
    ! of sin(theta) instead
    3 g1 = -(ct * mp1 * g2 + ss * t * g3)/d
    f1 = (mp1 * (2.0d0 * s * g2-ct * f2)-t * ss * (f3 + cot * g3))/d - cot * g1
    x(mp1) = g1
    dx(mp1) = f1
    if (m == 0) goto 4
    mp1 = m
    m = m - 1
    v = v + 1.0d0
    y = y - 1.0d0
    t = d
    d = dsqrt(v * y)
    g3 = g2
    g2 = g1
    f3 = f2
    f2 = f1
    goto 3
    ! explicit conversion to integer added
    4 maxsin = int(-72.0d0 / log10(s))
    
    ! maxsin is the max exponent of sin(theta) without underflow
    lpsafe = min0(lp1, maxsin)
    stom = 1.0d0
    fac = sign(1.0d0, dble((l / 2) * 2 - l)  +  0.50d0)
    
    ! multiply xlm by sin *  * m
    do m = 1, lpsafe
        x(m) = fac * x(m) * stom
        dx(m) = fac * dx(m) * stom
        stom = stom * s
    enddo
    
    ! set any remaining xlm to zero
    if (maxsin <= l) then
        mmm = maxsin + 1
        do m = mmm, lp1
            x(m) = 0.0d0
            dx(m) = 0.0d0
        enddo
    endif
    
    s = sos
    l = lsave
    
end subroutine lgndr


end module s20rts
