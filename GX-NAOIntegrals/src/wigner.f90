module wigner

    ! ============================================================================
    ! MIT License
    !
    ! Copyright (c) 2022 Oliver C. Gorton (ogorton@sdsu.edu)
    ! 
    ! Permission is hereby granted, free of charge, to any person obtaining a copy
    ! of this software and associated documentation files (the "Software"), to deal
    ! in the Software without restriction, including without limitation the rights
    ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    ! copies of the Software, and to permit persons to whom the Software is
    ! furnished to do so, subject to the following conditions:
    ! 
    ! The above copyright notice and this permission notice shall be included in all
    ! copies or substantial portions of the Software.
    ! 
    ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    ! SOFTWARE.
    ! ============================================================================

    ! Oliver Gorton, San Diego State University, 2021.12.7
    !
    ! Library of functions for computation of Wigner 3-j, 6-j and 9-j symbols
    ! using algebraic expressions in terms of factorials. Should be accurate
    ! to 10^{-10} relative error for values less than about j=20.
    !
    ! For an analysis of relative error compared to more modern methods, see
    ! arXiv:1504.08329 by H. T. Johansson and C. Forssen. A more accurate but
    ! slower method involves prime factorization of integers. In old Fortran,
    ! see work by Liqiang Wei:
    ! Computer Physics Communications 120 (1999) 222-230.
    !
    ! All integer arguments are 2*j in order to accomadate half-integer
    ! arguments while taking advantage of faster integer-arithmetic.
    ! Invalid arguments return 0d0 and program continues.
    !
    ! Optionally, compile with OpenMP to accelerate table initialization.
    !
    ! List of real(kind=8) functions:
    !     logfac(n)
    !     logdoublefac(n)
    !     triangle(two_j1, two_j2, two_j3)
    !     vector_couple(two_j1, two_m1, two_j2, two_m2, two_jc, two_mc)
    !     threej(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)
    !     threej_lookup(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
    !     sixj(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
    !     sixj_lookup(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
    !     ninej(two_j1,two_j2,two_j3,two_j4,two_j5,two_j6,two_j7,two_j8,two_j9)
    !
    ! List of subroutines:
    !     threej_table_init(min2j, max2j)
    !     sixj_table_init(min2j, max2j)

    implicit none
    real(kind=8), allocatable :: threej_table(:,:,:,:,:,:)
    real(kind=8), allocatable :: sixj_table(:,:,:,:,:,:)

    ! Default min/max values of 3-j, 6-j (x2) table arguments to store in
    ! lookup table when calling _table_init().
    integer :: tablemin2j = 0
    integer :: tablemax2j = 12

    ! These variables are modified by the lookup table functions to track
    ! the largest 2j values used during the program. Use these to inform future
    ! table sizes.
    integer :: tablemin_used
    integer :: tablemax_used

    contains

        function logfac(n)

            ! Computes log(n!)
            ! Log-factorial is used to delay numerical overflow.

            implicit none
            integer :: n
            real(kind=8) :: logfac
            integer :: i
            if (n<=0) then
                logfac = 0d0
                return
            end if
            logfac = 1d0
            do i = 1, n
                logfac = logfac * dble(i)
            end do
            logfac = log(logfac)
            return

        end function logfac

        function logdoublefac(n)

            ! Computes log(n!!)
            ! The double-factorial is NOT (n!)!. The double factorial is
            ! defined as the product of all the integers from 1 to n that
            ! have the same parity (even or odd) as n.

            implicit none
            integer :: n
            real(kind=8) :: logdoublefac
            integer :: i, imin

            logdoublefac = 0d0

            if (n <= 0) return

            if (mod(n,2) == 0) then
                ! n is even
                imin = 2
            else
                ! n is odd
                imin = 1
            end if

            logdoublefac = 1d0
            do i = imin, n, 2
                logdoublefac = logdoublefac * dble(i)
            end do

            logdoublefac = log(logdoublefac)

            return

        end function logdoublefac

        function triangle(two_j1, two_j2, two_j3) result(delta)

            ! Computes the triangle functions, typically denoted as
            !     \delta(abc) = (a+b-c)!(a-b+c)!(b+c-a)!/(a+b+c+1)!

            implicit none
            integer :: two_j1, two_j2, two_j3
            integer :: c1, c2, c3, c4
            real(kind=8) :: delta

            delta = 0.0D0
            c1 = two_j1 + two_j2 - two_j3
            c2 = two_j1 - two_j2 + two_j3
            c3 =-two_j1 + two_j2 + two_j3
            c4 = two_j1 + two_j2 + two_j3

            if (c1<0) return
            if (c2<0) return
            if (c3<0) return
            if (c4+1<0) return

            if (mod(c1,2)/=0) return
            if (mod(c2,2)/=0) return
            if (mod(c3,2)/=0) return
            if (mod(c4,2)/=0) return

            delta = exp(0.5d0*(logfac(c1/2)+logfac(c2/2)&
                +logfac(c3/2)-logfac(c4/2+1)))

            return

        end function triangle

        function vector_couple(two_j1, two_m1, two_j2, &
                               two_m2, two_jc, two_mc) result(cg)

            ! Computes the Clebsh-Gordon vector-coupling coefficient
            !  (j1 m1 j2 m2 | j1 j2 jc mc)
            ! using algebraic expressions in factorials from Edmonds.

            implicit none
            integer :: two_j1,two_j2,two_jc,two_m1,two_m2,two_mc
            real(kind=8) :: cg, fac_prod, fac_sum, den
            integer :: t1, t2, t3, t4
            integer :: d1, d2, d3, d4, d5, d6
            integer :: dd1, dd2
            integer :: zmin, zmax, z

            cg = 0d0
            if (two_m1+two_m2 /= two_mc) return
            if (two_j1<0) return
            if (two_j2<0) return
            if (two_jc<0) return
            if (abs(two_m1)>two_j1) return
            if (abs(two_m2)>two_j2) return
            if (abs(two_jc)>two_jc) return
            if (mod(two_j1+two_m1,2) /= 0) return
            if (mod(two_j2+two_m2,2) /= 0) return
            if (mod(two_jc+two_mc,2) /= 0) return

            t1 = ( two_j1 + two_j2 - two_jc)/2
            t2 = ( two_j1 - two_j2 + two_jc)/2
            t3 = (-two_j1 + two_j2 + two_jc)/2
            t4 = ( two_j1 + two_j2 + two_jc)/2
            if (t1<0) return
            if (t2<0) return
            if (t3<0) return

            d1 = (two_j1 + two_m1)/2
            d2 = (two_j1 - two_m1)/2
            d3 = (two_j2 + two_m2)/2
            d4 = (two_j2 - two_m2)/2
            d5 = (two_jc + two_mc)/2
            d6 = (two_jc - two_mc)/2

            dd1 = (two_jc - two_j2 + two_m1)/2
            dd2 = (two_jc - two_j1 - two_m2)/2

            fac_prod = sqrt(dble(two_jc)+1d0) &
                * triangle(two_j1,two_j2,two_jc) &
                * exp(0.5d0 * (logfac(d1) + logfac(d2) &
                              +logfac(d3) + logfac(d4) &
                              +logfac(d5) + logfac(d6)))

            zmin = max(0, -dd1, -dd2)
            zmax = min(t1, d2, d3)

            fac_sum = 0d0
            do z = zmin, zmax
                den =- (logfac(z) + logfac(t1-z) + logfac(d2-z) &
                      + logfac(d3-z) + logfac(dd1+z) + logfac(dd2+z))
                fac_sum = fac_sum + (-1) ** (z) * exp(den)
            end do

            cg = fac_prod * fac_sum

            return

        end function vector_couple

        function threej(two_j1, two_j2, two_j3,&
                        two_m1, two_m2, two_m3) result(tj)

            ! Computes the Wigner 3-J symbol with arguments
            !   two_j1/2 two_j1/2 two_j3/2
            !   two_m1/2 two_m2/2 two_m3/2
            ! using clebsh-gordon vector-coupling coefficients.

            implicit none
            integer :: two_j1,two_j2,two_j3,two_m1,two_m2,two_m3
            real(kind=8) :: tj

            tj = (-1) **((two_j1 - two_j2 - two_m3)/2)/sqrt(dble(two_j3)+1d0) &
                * vector_couple(two_j1, two_m1, two_j2, two_m2, two_j3, -two_m3)

        end function threej

        subroutine threej_table_init(min2j, max2j)

            ! Initializes the table of 3-j symbols into memory accessible from
            ! the module. If the optional arguments are given, they override the
            ! default values set in the module variables for the min/max j-
            ! values stored in the table.

            implicit none
            integer, optional :: min2j, max2j
            integer :: a, b
            integer :: i, j, k, l, m, n

            integer(kind=8) :: ti, tf, clock_rate
            real(kind=8) :: dummy_real

            call system_clock(count_rate = clock_rate)
            call system_clock(count = ti)

            print '(a)',"Initializing three-j symbol table..."
            if (.not. present(min2j)) then
                a = tablemin2j
            else
                a = min2j
                tablemin2j = a
            end if
            if (.not. present(max2j)) then
                b = tablemax2j
            else
                b = max2j
                tablemax2j = b
            end if
            print*,"Table min. 2J:",a
            print*,"Table max. 2J:",b
            print "(a,f10.2)",'Memory required (MB):',real(sizeof(dummy_real) &
                                    * (b - a + 1)**6d0) * 10d0 ** (-6)

            if (allocated(threej_table)) deallocate(threej_table)
            allocate(threej_table(a:b,a:b,a:b,a:b,a:b,a:b))
            threej_table = 0d0
            do n = a, b
              do m = a, b
                do l = a, b
                  do k = a, b
                    do j = a, b
                      do i = a, b
                        threej_table(i,j,k,l,m,n) = threej(i,j,k,l,m,n)
                      end do
                    end do
                  end do
                end do
              end do
            end do

            call system_clock(count = tf)
            print*,'Table has been saved to memory.'
            print*,'Seconds to initialize:',real((tf-ti))/real(clock_rate)

        end subroutine threej_table_init

        function threej_lookup(two_j1, two_j2, two_j3,&
                               two_l1, two_l2, two_l3) result(tj)

            ! Function to lookup a value of the 3-j symbol stored in the
            ! lookup table. If the requested value has any argument outside
            ! the defined lookup table values, it is computed using the 3-j
            ! function instead.
            ! This function also updates module variables tracking the largest
            ! j values requested from the lookup function, to inform future
            ! table initializations.


            implicit none
            integer :: two_j1,two_j2,two_j3,two_l1,two_l2,two_l3
            real(kind=8) :: tj
            integer :: minj, maxj

            minj = min(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            maxj = max(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            tablemin_used = min(tablemin_used, minj)
            tablemax_used = max(tablemax_used, maxj)

            if (minj < tablemin2j .or. maxj > tablemax2j) then
                tj = threej(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            else
                tj = threej_table(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            end if

            return

        end function threej_lookup

        function sixj(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3) result(sj)

            ! Computes the wigner six-j symbol with arguments
            !    two_j1/2 two_j2/2 two_j3/2
            !    two_l1/2 two_l2/2 two_l3/2
            ! using explicit algebraic expressions from Edmonds (1955/7).

            implicit none
            integer :: two_j1,two_j2,two_j3,two_l1,two_l2,two_l3
            integer :: z, zmin, zmax
            integer  :: n2, n3, n4
            integer  :: d1, d2, d3, d4
            real(kind=8) :: sj
            real(kind=8) :: fac_sum, triangle_prod, num, den

            sj = 0d0
            if ( two_j1 < 0) return
            if ( two_j2 < 0) return
            if ( two_j3 < 0) return
            if ( two_l1 < 0) return
            if ( two_l2 < 0) return
            if ( two_l3 < 0) return
            if ( two_j1 < abs(two_j2-two_j3)) return
            if ( two_j1 > two_j2+two_j3 ) return
            if ( two_j1 < abs(two_l2-two_l3)) return
            if ( two_j1 > two_l2+two_l3 ) return
            if ( two_l1 < abs(two_j2-two_l3)) return
            if ( two_l1 > two_j2+two_l3 ) return
            if ( two_l1 < abs(two_l2-two_j3)) return
            if ( two_l1 > two_l2+two_j3 ) return

            n2 = (two_j1+two_j2+two_l1+two_l2)/2
            n3 = (two_j2+two_j3+two_l2+two_l3)/2
            n4 = (two_j3+two_j1+two_l3+two_l1)/2

            d1 = (two_j1+two_j2+two_j3)/2
            d2 = (two_j1+two_l2+two_l3)/2
            d3 = (two_l1+two_j2+two_l3)/2
            d4 = (two_l1+two_l2+two_j3)/2

            triangle_prod = triangle(two_j1,two_j2,two_j3) &
                           *triangle(two_j1,two_l2,two_l3) &
                           *triangle(two_l1,two_j2,two_l3) &
                           *triangle(two_l1,two_l2,two_j3)

            zmin = max(d1, d2, d3, d4)
            zmax = min(n2, n3, n4)

            fac_sum = 0d0

            do z = zmin, zmax
                num = logfac(z+1)
                den = logfac(n2-z) + logfac(n3-z) + logfac(n4-z) &
                    + logfac(z-d1) + logfac(z-d2) + logfac(z-d3) + logfac(z-d4)
                fac_sum = fac_sum + (-1) ** (z) * exp(num - den)
            end do

            sj = triangle_prod * fac_sum

           return

        end function sixj

        subroutine sixj_table_init(min2j, max2j)

            implicit none
            integer, optional :: min2j, max2j
            integer :: a, b
            integer :: i, j, k, l, m, n

            integer(kind=8) :: ti, tf, clock_rate
            real(kind=8) :: dummy_real

            call system_clock(count_rate = clock_rate)
            call system_clock(count = ti)

            print '(a)',"Initializing six-j symbol table..."
            if (.not. present(min2j)) then
                a = tablemin2j
            else
                a = min2j
                tablemin2j = a
            end if
            if (.not. present(max2j)) then
                b = tablemax2j
            else
                b = max2j
                tablemax2j = b
            end if
            print*,"Table min. 2J:",a
            print*,"Table max. 2J:",b
            print "(a,f10.2)","Memory required (MB):",real(sizeof(dummy_real) &
                                    * (b - a + 1)**6d0,4)*10d0**(-6)

            if (allocated(sixj_table)) deallocate(sixj_table)
            allocate(sixj_table(a:b,a:b,a:b,a:b,a:b,a:b))
            sixj_table = 0d0

            do n = a, b
              do m = a, b
                do l = a, b
                  do k = a, b
                    do j = a, b
                      do i = a, b
                        sixj_table(i,j,k,l,m,n) = sixj(i,j,k,l,m,n)
                      end do
                    end do
                  end do
                end do
              end do
            end do

            call system_clock(count = tf)
            print '(a)', 'Table has been saved to memory.'
            print '(a, f10.4)', 'Seconds to initialize:', real((tf-ti))/real(clock_rate)

        end subroutine sixj_table_init

        function sixj_lookup(two_j1, two_j2, two_j3,&
                             two_l1, two_l2, two_l3) result(sj)

            implicit none
            integer :: two_j1,two_j2,two_j3,two_l1,two_l2,two_l3
            real(kind=8) :: sj
            integer :: minj, maxj

            minj = min(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            maxj = max(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            tablemin_used = min(tablemin_used, minj)
            tablemax_used = max(tablemax_used, maxj)

            if (minj < tablemin2j .or. maxj > tablemax2j) then
                sj = sixj(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            else
                sj = sixj_table(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3)
            end if

            return

        end function sixj_lookup

        function ninej(two_j1, two_j2, two_j3,&
                       two_j4, two_j5, two_j6,&
                       two_j7, two_j8, two_j9) result(nj)

            implicit none
            integer :: two_j1,two_j2,two_j3
            integer :: two_j4,two_j5,two_j6
            integer :: two_j7,two_j8,two_j9
            integer :: z, zmin, zmax
            real(kind=8) :: nj

            nj = 0d0

            if (two_j1 < 0) return
            if (two_j2 < 0) return
            if (two_j3 < 0) return
            if (two_j4 < 0) return
            if (two_j5 < 0) return
            if (two_j6 < 0) return
            if (two_j7 < 0) return
            if (two_j8 < 0) return
            if (two_j9 < 0) return

            zmin = max(abs(two_j1-two_j9),abs(two_j4-two_j8),abs(two_j2-two_j6))
            zmax = min(two_j1+two_j9, two_j4+two_j8, two_j2+two_j6)

            do z = zmin, zmax
                nj = nj + (-1) ** z * (z + 1d0) &
                    * sixj_lookup(two_j1,two_j4,two_j7,two_j8,two_j9,     z) &
                    * sixj_lookup(two_j2,two_j5,two_j8,two_j4,     z,two_j6) &
                    * sixj_lookup(two_j3,two_j6,two_j9,     z,two_j1,two_j2)
            end do

            return

        end function ninej

end module wigner