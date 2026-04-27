MODULE mod_gn

implicit none

CONTAINS

      SUBROUTINE gn(iplc, xi, N)

! FORTRAN version of gn.m

!  subroutine gn(iplc, xi, N)

!  Computation of shape functions.

!  input:
!    iplc     :  interpolation code
!    xi       :  local coordinates of integration point

!  output:
!    N        :  array of shape functions


!  Variable declaration
      implicit none

! INPUT --------------------------------------------------------------------------------
      character(len=3),intent(in)                                  :: iplc
      real*8,intent(in)                                            :: xi(*)

! OUTPUT -------------------------------------------------------------------------------
      real*8,intent(out)                                           :: N(*)

!  local coordinates -------------------------------------------------------------------
      real*8                                                        :: xi1, xi2, xi3, xi4

!  Working variables -------------------------------------------------------------------
      character(len=80)                                             :: errmsgtxt
      real*8                                                        :: dw1, dw2, dw3, dw4, dw5, dw6


!   External routines

!   MATLAB matrix manipulation
 !     external          mexerrmsgtxt


!   Executable statements


!   One-dimensional
!   ---------------

      if ((iplc(:1) .eq. 'b').OR.(iplc(:1) .eq. 'l')) then

!       local coordinates of integration point
        xi1 = xi(1)

!       2 nodes
        if ((iplc .eq. 'b02').OR.(iplc .eq. 'l02')) then
          N(1) =  0.5d0 * (1.0d0 - xi1)
          N(2) =  0.5d0 * (1.0d0 + xi1)

!       3 nodes
        elseif ((iplc .eq. 'b03').OR.(iplc .eq. 'l03')) then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          N(1) = -0.5d0 * xi1 * dw2
          N(2) =          dw1 * dw2
          N(3) =  0.5d0 * dw1 * xi1

!       4 nodes
        elseif ((iplc .eq. 'b03').OR.(iplc .eq. 'l03')) then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          N(1)= 0.5d0 * dw2 - 0.5d0 * (1.0d0 - xi1**2.0d0) + (-9.0d0 * xi1**3.0d0 + xi1**2.0d0 + 9.0d0 * xi1 - 1.0d0)/16.0d0
          N(2)= (1.0d0 - xi1*xi1) + (27.0d0 * xi1**3.0d0 + 7.0d0 * xi1*xi1 - 27.0d0 * xi1 - 7.0d0)/16.0d0
          N(3)= (-27.0d0 * xi1**3.0d0 - 9.0d0 * xi1 * xi1 + 27.0d0 * xi1 + 9.0d0)/16.0d0
          N(4)= 0.5*dw1 - 0.5d0 * (1.0d0 - xi1**2.0d0) + (9.0d0 * xi1**3.0d0 + xi1**2.0d0 - 9.0d0 * xi1 - 1.0d0)/16.0d0

!       illegal interpolation code
        else
          write(errmsgtxt, 9960) iplc
 9960     format('Illegal interpolation code: ', a3, '.')
        !  call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif


!  Two-dimensional, triangular
!  ---------------------------

      elseif (iplc(:1) .eq. 't') then

!       local coordinates of integration point
        xi1 = xi(1)
        xi2 = xi(2)
        xi3 = xi(3)

!       3 nodes
        if (iplc .eq. 't03') then
          N(1) =  xi1
          N(2) =  xi2
          N(3) =  xi3

!       6 nodes
        elseif (iplc .eq. 't06') then
          N(1) =          xi1 * (2.0d0 * xi1 - 1.0d0)
          N(2) =  4.0d0 * xi1 *    xi2
          N(3) =          xi2 * (2.0d0 * xi2 - 1.0d0)
          N(4) =  4.0d0 * xi2 *    xi3
          N(5) =          xi3 * (2.0d0 * xi3 - 1.0d0)
          N(6) =  4.0d0 * xi1 *    xi3

!       illegal interpolation code
        else
          write(errmsgtxt, 9970) iplc
 9970     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif


!  Two-dimensional, quadrilateral
!  ------------------------------

      elseif (iplc(:1) .eq. 'q') then

!       local coordinates of integration point
        xi1 = xi(1)
        xi2 = xi(2)

!       4 nodes
        if (iplc .eq. 'q04') then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          N(1) =  0.25d0 * dw2 * dw4
          N(2) =  0.25d0 * dw1 * dw4
          N(3) =  0.25d0 * dw1 * dw3
          N(4) =  0.25d0 * dw2 * dw3

!       8 nodes
        elseif (iplc .eq. 'q08') then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          N(1) = -0.25d0 * dw2 * dw4 * (1+xi1+xi2)
          N(2) =  0.50d0 * dw1 * dw2 * dw4
          N(3) = -0.25d0 * dw1 * dw4 * (1-xi1+xi2)
          N(4) =  0.50d0 * dw1 * dw3 * dw4
          N(5) = -0.25d0 * dw1 * dw3 * (1-xi1-xi2)
          N(6) =  0.50d0 * dw1 * dw2 * dw3
          N(7) = -0.25d0 * dw2 * dw3 * (1+xi1-xi2)
          N(8) =  0.50d0 * dw2 * dw3 * dw4

!       9 nodes
        elseif (iplc .eq. 'q09') then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          N(1) =  0.25d0 * xi1 * dw2 * xi2 * dw4
          N(2) = -0.50d0 * dw1 * dw2 * xi2 * dw4
          N(3) = -0.25d0 * dw1 * xi1 * xi2 * dw4
          N(4) =  0.50d0 * dw1 * xi1 * dw3 * dw4
          N(5) =  0.25d0 * dw1 * xi1 * dw3 * xi2
          N(6) =  0.50d0 * dw1 * dw2 * dw3 * xi2
          N(7) = -0.25d0 * xi1 * dw2 * dw3 * xi2
          N(8) = -0.50d0 * xi1 * dw2 * dw3 * dw4
          N(9) =           dw1 * dw2 * dw3 * dw4

!       12 nodes
        elseif (iplc .eq. 'q12') then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          N(1) = 1.0d0/32.0d0 * dw2 * dw4 * (-10.0d0 + 9.0d0 * (xi1 * xi1 + xi2 * xi2))
          N(2) = 9.0d0/32.0d0 * dw4 * (1.0d0 -xi1 * xi1) * (1.0d0 - 3.0d0 * xi1)
          N(3) = 9.0d0/32.0d0 * dw4 * (1.0d0 -xi1 * xi1) * (1.0d0 + 3.0d0 * xi1)
          N(4) = 1.0d0/32.0d0 * dw1 * dw4 * (-10.0d0 + 9.0d0 * (xi1 * xi1 + xi2 * xi2))
          N(5) = 9.0d0/32.0d0 * dw1 * (1.0d0 -xi2 * xi2) * (1.0d0 - 3.0d0 * xi2)
          N(6) = 9.0d0/32.0d0 * dw1 * (1.0d0 -xi2 * xi2) * (1.0d0 + 3.0d0 * xi2)
          N(7) = 1.0d0/32.0d0 * dw1 * dw3 * (-10.0d0 + 9.0d0 * (xi1 * xi1 + xi2 * xi2))
          N(8) = 9.0d0/32.0d0 * dw3 * (1.0d0 -xi1 * xi1) * (1.0d0 + 3.0d0 * xi1)
          N(9) = 9.0d0/32.0d0 * dw3 * (1.0d0 -xi1 * xi1) * (1.0d0 - 3.0d0 * xi1)
          N(10)= 1.0d0/32.0d0 * dw2 * dw3 * (-10.0d0 + 9.0d0 * (xi1 * xi1 + xi2 * xi2))
          N(11)= 9.0d0/32.0d0 * dw2 * (1.0d0 -xi2 * xi2) * (1.0d0 + 3.0d0 * xi2)
          N(12)= 9.0d0/32.0d0 * dw2 * (1.0d0 -xi2 * xi2) * (1.0d0 - 3.0d0 * xi2)

!       illegal interpolation code
        else
          write(errmsgtxt, 9980) iplc
 9980     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif


!  Three-dimensional, quadrilateral
!  --------------------------------

      elseif (iplc(:1) .eq. 'v') then

!       local coordinates of integration point
        xi1 = xi(1)
        xi2 = xi(2)
        xi3 = xi(3)

!       8 nodes
        if (iplc .eq. 'v08') then
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          dw5 = 1.0d0 + xi3
          dw6 = 1.0d0 - xi3
          N(1) =  0.125d0 * dw2 * dw4 * dw6
          N(2) =  0.125d0 * dw1 * dw4 * dw6
          N(3) =  0.125d0 * dw1 * dw3 * dw6
          N(4) =  0.125d0 * dw2 * dw3 * dw6
          N(5) =  0.125d0 * dw2 * dw4 * dw5
          N(6) =  0.125d0 * dw1 * dw4 * dw5
          N(7) =  0.125d0 * dw1 * dw3 * dw5
          N(8) =  0.125d0 * dw2 * dw3 * dw5

!       20 nodes
        elseif (iplc .eq. 'v20') then
        N(1)  = 0.125d0 * (1.0d0-xi1) * (1.0d0-xi2) * (1.0d0-xi3)* (-2.0d0-xi1-xi2-xi3)
        N(2)  = 0.25d0  * (1.0d0-xi2) * (1.0d0-xi3) * (1.0d0-xi1**2.0d0)
        N(3)  = 0.125d0 * (1.0d0+xi1) * (1.0d0-xi2) * (1.0d0-xi3)* (-2.0d0+xi1-xi2-xi3)
        N(4)  = 0.25d0  * (1.0d0+xi1) * (1.0d0-xi3) * (1.0d0-xi2**2.0d0)
        N(5)  = 0.125d0 * (1.0d0+xi1) * (1.0d0+xi2) * (1.0d0-xi3)* (-2.0d0+xi1+xi2-xi3)
        N(6)  = 0.25d0  * (1.0d0+xi2) * (1.0d0-xi3) * (1.0d0-xi1**2.0d0)
        N(7)  = 0.125d0 * (1.0d0-xi1) * (1.0d0+xi2) * (1.0d0-xi3)* (-2.0d0-xi1+xi2-xi3)
        N(8)  = 0.25d0  * (1.0d0-xi1) * (1.0d0-xi3) * (1.0d0-xi2**2.0d0)
        N(9)  = 0.25d0  * (1.0d0-xi1) * (1.0d0-xi2) * (1.0d0-xi3**2.0d0)
        N(10) = 0.25d0  * (1.0d0+xi1) * (1.0d0-xi2) * (1.0d0-xi3**2.0d0)
        N(11) = 0.25d0  * (1.0d0+xi1) * (1.0d0+xi2) * (1.0d0-xi3**2.0d0)
        N(12) = 0.25d0  * (1.0d0-xi1) * (1.0d0+xi2) * (1.0d0-xi3**2.0d0)
        N(13) = 0.125d0 * (1.0d0-xi1) * (1.0d0-xi2) * (1.0d0+xi3)* (-2.0d0-xi1-xi2+xi3)
        N(14) = 0.25d0  * (1.0d0-xi2) * (1.0d0+xi3) * (1.0d0-xi1**2.0d0)
        N(15) = 0.125d0 * (1.0d0+xi1) * (1.0d0-xi2) * (1.0d0+xi3)* (-2.0d0+xi1-xi2+xi3)
        N(16) = 0.25d0  * (1.0d0+xi1) * (1.0d0+xi3) * (1.0d0-xi2**2.0d0)
        N(17) = 0.125d0 * (1.0d0+xi1) * (1.0d0+xi2) * (1.0d0+xi3)* (-2.0d0+xi1+xi2+xi3)
        N(18) = 0.25d0  * (1.0d0+xi2) * (1.0d0+xi3) * (1.0d0-xi1**2.0d0)
        N(19) = 0.125d0 * (1.0d0-xi1) * (1.0d0+xi2) * (1.0d0+xi3)* (-2.0d0-xi1+xi2+xi3)
        N(20) = 0.25d0  * (1.0d0-xi1) * (1.0d0+xi3) * (1.0d0-xi2**2.0d0)

!       illegal interpolation code
        else
          write(errmsgtxt, 9985) iplc
 9985     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif

!  Three-dimensional, tetrahedral
!  --------------------------------

      elseif (iplc(:1) .eq. 'h') then

!       local coordinates of integration point
        xi1 = xi(1)
        xi2 = xi(2)
        xi3 = xi(3)
        xi4 = xi(4)

!       4 nodes
        if (iplc .eq. 'h04') then
          N(1) =  xi1
          N(2) =  xi2
          N(3) =  xi3
          N(4) =  xi4

!       10 nodes
        elseif (iplc .eq. 'h10') then
        N(1) = (2.0d0 * xi1 - 1.0d0) * xi1
        N(2) = (2.0d0 * xi2 - 1.0d0) * xi2
        N(3) = (2.0d0 * xi3 - 1.0d0) * xi3
        N(4) = (2.0d0 * xi4 - 1.0d0) * xi4
        N(5) = 4.0d0 * xi1 * xi2
        N(6) = 4.0d0 * xi2 * xi3
        N(7) = 4.0d0 * xi1 * xi3
        N(8) = 4.0d0 * xi1 * xi4
        N(9) = 4.0d0 * xi2 * xi4
        N(10)= 4.0d0 * xi3 * xi4

!       illegal interpolation code
        else
          write(errmsgtxt, 9990) iplc
 9990     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif


! Illegal interpolation code
! --------------------------

      else
        write(errmsgtxt, 9995) iplc
 9995   format('Illegal interpolation code: ', a3, '.')
!        call mexerrmsgtxt(errmsgtxt)
        write(*,*) "Illegal interpolation code"
        STOP
      endif

       END SUBROUTINE gn
END MODULE mod_gn
