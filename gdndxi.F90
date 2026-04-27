MODULE mod_gdndxi

implicit none

CONTAINS

      SUBROUTINE gdndxi(iplc, xi, nnode, dNdxi)

! FORTRAN version of gdndxi.m

!  subroutine gdndxi(iplc, xi, dNdxi, ld)

!  Calculation of derivatives of shape functions with respect to
!  local coordinates.

!  input:
!    iplc     :  interpolation code
!    xi       :  local coordinates of integration point
!    ld       :  leading dimension of dNdxi-array

!  output:
!    dNdxi    :  array of derivatives

!   Variable declaration
      implicit none

! INPUT -------------------------------------------------------------------------
      character(len=3),intent(in)                               :: iplc
      integer*8,intent(in)                                      :: nnode
      real*8,intent(in)                                         :: xi(*)

! OUTPUT ------------------------------------------------------------------------
      real*8,intent(out)                                        :: dNdxi(nnode, *)

!  Local coordinates ------------------------------------------------------------
      real*8                                                    :: xi1, xi2, xi3, xi4

!  Working variables ------------------------------------------------------------
      character(len=80)                                         :: errmsgtxt
      real*8                                                    :: dw1, dw2, dw3, dw4, dw5, dw6
 character(len=10) :: test 

! External routines

!  MATLAB matrix manipulation
!     external          mexerrmsgtxt


! Executable statements


! One-dimensional
! ---------------

      if ((iplc(:1) .eq. 'b').OR.(iplc(:1) .eq. 'l')) then

!       local coordinates of integration point
        xi1 = xi(1)

!       2 nodes
        if ((iplc .eq. 'b02').OR.(iplc .eq. 'l02')) then
          dNdxi(1, 1) = -0.5d0
          dNdxi(2, 1) =  0.5d0

!       3 nodes
        elseif ((iplc .eq. 'b03').OR.(iplc .eq. 'l03')) then
          dNdxi(1, 1) =  xi1 - 0.5d0
          dNdxi(2, 1) = -2.0d0 * xi1
          dNdxi(3, 1) =  xi1 + 0.5d0

!       4 nodes
        elseif ((iplc .eq. 'b04').OR.(iplc .eq. 'l04')) then
          dNdxi(1, 1) = -0.5d0 + xi1 + (-27.0d0 * xi1 * xi1 + 2 * xi1 + 9.0d0)/16.0d0
          dNdxi(2, 1) = -2.0d0 * xi1 + (81.0d0 * xi1 * xi1 + 14.0d0 * xi1 - 27.0d0)/16.0d0
          dNdxi(3, 1) = (-81.0d0 * xi1 * xi1 - 18.0d0 * xi1 - 27.0d0)/16.0d0
          dNdxi(4, 1) = 0.5d0 + xi1 + (27.0d0 * xi1 * xi1 + 2 * xi1 - 9.0d0)/16.0d0

!       illegal interpolation code
        else
          write(errmsgtxt, 9960) iplc
 9960     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
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
          dNdxi(1, 1) =  1.0d0
          dNdxi(2, 1) =  0.0d0
          dNdxi(3, 1) = -1.0d0
          dNdxi(1, 2) =  0.0d0
          dNdxi(2, 2) =  1.0d0
          dNdxi(3, 2) = -1.0d0

!       6 nodes
        elseif (iplc .eq. 't06') then

          dNdxi(1, 1:2) =[  4.0d0 * xi1 - 1.0d0      ,   0.0d0                     ]
          dNdxi(2, 1:2) =[  4.0d0 * xi2              ,   4.0d0 * xi1               ]
          dNdxi(3, 1:2) =[  0.0d0                    ,   4.0d0 * xi2 - 1.0d0       ]
          dNdxi(4, 1:2) =[ -4.0d0 * xi2              ,   4.0d0 * xi3 - 4.0d0 * xi2 ]
          dNdxi(5, 1:2) =[ -4.0d0 * xi3 + 1.0d0      ,   -4.0d0 * xi3 + 1.0d0      ]
          dNdxi(6, 1:2) =[ 4.0d0 * xi3 - 4.0d0 * xi1,   -4.0d0 * xi1               ]
!
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
          dw1 = 0.25d0 * (1.0d0 + xi1)
          dw2 = 0.25d0 * (1.0d0 - xi1)
          dw3 = 0.25d0 * (1.0d0 + xi2)
          dw4 = 0.25d0 * (1.0d0 - xi2)
          dNdxi(1, 1) = -dw4
          dNdxi(2, 1) =  dw4
          dNdxi(3, 1) =  dw3
          dNdxi(4, 1) = -dw3
          dNdxi(1, 2) = -dw2
          dNdxi(2, 2) = -dw1
          dNdxi(3, 2) =  dw1
          dNdxi(4, 2) =  dw2

!       8 nodes
        elseif (iplc .eq. 'q08') then
          dw1 = 1.0d0 + xi(1)
          dw2 = 1.0d0 - xi(1)
          dw3 = 1.0d0 + xi(2)
          dw4 = 1.0d0 - xi(2)
          dNdxi(1, 1) =  0.25d0 * dw4 * (2.0d0 * xi1 + xi2)
          dNdxi(2, 1) =          -dw4 * xi1
          dNdxi(3, 1) =  0.25d0 * dw4 * (2.0d0 * xi1 - xi2)
          dNdxi(4, 1) =  0.50d0 * dw3 * dw4
          dNdxi(5, 1) =  0.25d0 * dw3 * (2.0d0 * xi1 + xi2)
          dNdxi(6, 1) =          -dw3 * xi1
          dNdxi(7, 1) =  0.25d0 * dw3 * (2.0d0 * xi1 - xi2)
          dNdxi(8, 1) = -0.50d0 * dw3 * dw4
          dNdxi(1, 2) =  0.25d0 * dw2 * (2.0d0 * xi2 + xi1)
          dNdxi(2, 2) = -0.50d0 * dw1 * dw2
          dNdxi(3, 2) =  0.25d0 * dw1 * (2.0d0 * xi2 - xi1)
          dNdxi(4, 2) =          -dw1 * xi2
          dNdxi(5, 2) =  0.25d0 * dw1 * (2.0d0 * xi2 + xi1)
          dNdxi(6, 2) =  0.50d0 * dw1 * dw2
          dNdxi(7, 2) =  0.25d0 * dw2 * (2.0d0 * xi2 - xi1)
          dNdxi(8, 2) =          -dw2 * xi2

!       9 nodes
        elseif (iplc .eq. 'q09') then
          dw1 = 1.0d0 + 2.0d0 * xi1
          dw2 = 1.0d0 - 2.0d0 * xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          dNdxi(1, 1) =  0.25d0 * dw2 * xi2 * dw4
          dNdxi(2, 1) =           xi1 * xi2 * dw4
          dNdxi(3, 1) = -0.25d0 * dw1 * xi2 * dw4
          dNdxi(4, 1) =  0.50d0 * dw1 * dw3 * dw4
          dNdxi(5, 1) =  0.25d0 * dw1 * dw3 * xi2
          dNdxi(6, 1) =          -xi1 * dw3 * xi2
          dNdxi(7, 1) = -0.25d0 * dw2 * dw3 * xi2
          dNdxi(8, 1) = -0.50d0 * dw2 * dw3 * dw4
          dNdxi(9, 1) = -2.00d0 * xi1 * dw3 * dw4
          dw1 = 1.0d0 + xi1
          dw2 = 1.0d0 - xi1
          dw3 = 1.0d0 + 2.0d0 * xi2
          dw4 = 1.0d0 - 2.0d0 * xi2
          dNdxi(1, 2) =  0.25d0 * xi1 * dw2 * dw4
          dNdxi(2, 2) = -0.50d0 * dw1 * dw2 * dw4
          dNdxi(3, 2) = -0.25d0 * dw1 * xi1 * dw4
          dNdxi(4, 2) =          -dw1 * xi1 * xi2
          dNdxi(5, 2) =  0.25d0 * dw1 * xi1 * dw3
          dNdxi(6, 2) =  0.50d0 * dw1 * dw2 * dw3
          dNdxi(7, 2) = -0.25d0 * xi1 * dw2 * dw3
          dNdxi(8, 2) =           xi1 * dw2 * xi2
          dNdxi(9, 2) = -2.00d0 * dw1 * dw2 * xi2

!       12 nodes
        elseif (iplc .eq. 'q12') then
          dw1 = 1.0d0 + 3.0d0 * xi1
          dw2 = 1.0d0 - 3.0d0 * xi1
          dw3 = 1.0d0 + xi2
          dw4 = 1.0d0 - xi2
          dw5 = 1.0d0 + 3.0d0 * xi2
          dw6 = 1.0d0 - 3.0d0 * xi2
          dNdxi(1, 1) = 1.0d0/32.0d0 * dw4 * (-(-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi1 * (1 - xi1))
          dNdxi(2, 1) = 9.0d0/32.0d0 * dw4 * (-2.0d0 * xi1 * dw2 - 3.0d0 * (1- xi1 * xi1))
          dNdxi(3, 1) = 9.0d0/32.0d0 * dw4 * (-2.0d0 * xi1 * dw1 - 3.0d0 * (1- xi1 * xi1))
          dNdxi(4, 1) = 1.0d0/32.0d0 * dw4 * ( (-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi1 * (1 + xi1))
          dNdxi(5, 1) = 9.0d0/32.0d0 * 1.0d0 * (1-xi2*xi2) * dw6
          dNdxi(6, 1) = 9.0d0/32.0d0 * 1.0d0 * (1-xi2*xi2) * dw5
          dNdxi(7, 1) = 1.0d0/32.0d0 * dw3 * ( (-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi1 * (1 + xi1))
          dNdxi(8, 1) = 9.0d0/32.0d0 * dw3 * (-2.0d0 * xi1 * dw1 + 3.0d0 * (1- xi1 * xi1))
          dNdxi(9, 1) = 9.0d0/32.0d0 * dw3 * (-2.0d0 * xi1 * dw2 - 3.0d0 * (1- xi1 * xi1))
          dNdxi(10,1) = 1.0d0/32.0d0 * dw3 * (-(-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi1 * (1 - xi1))
          dNdxi(11,1) = 9.0d0/32.0d0 * (-1.0d0) * (1-xi2*xi2) * dw5
          dNdxi(12,1) = 9.0d0/32.0d0 * (-1.0d0) * (1-xi2*xi2) * dw6
          dw3 = 1.0d0 + xi1
          dw4 = 1.0d0 - xi1
          dNdxi(1, 2) = 1.0d0/32.0d0 * dw4 * (-(-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi2 * (1 - xi2))
          dNdxi(2, 2) = 9.0d0/32.0d0 * (-1.0d0) * (1-xi1*xi1) * dw2
          dNdxi(3, 2) = 9.0d0/32.0d0 * (-1.0d0) * (1-xi1*xi1) * dw1
          dNdxi(4, 2) = 1.0d0/32.0d0 * dw3 * (-(-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi2 * (1 - xi2))
          dNdxi(5, 2) = 9.0d0/32.0d0 * dw3 * (-2.0d0 * xi2 * dw6 + 3.0d0 * (1- xi2 * xi2))
          dNdxi(6, 2) = 9.0d0/32.0d0 * dw3 * (-2.0d0 * xi2 * dw5 + 3.0d0 * (1- xi2 * xi2))
          dNdxi(7, 2) = 1.0d0/32.0d0 * dw3 * ( (-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi2 * (1 + xi2))
          dNdxi(8, 2) = 9.0d0/32.0d0 * 1.0d0 * (1-xi1*xi1) * dw1
          dNdxi(9, 2) = 9.0d0/32.0d0 * 1.0d0 * (1-xi1*xi1) * dw2
          dNdxi(10,2) = 1.0d0/32.0d0 * dw4 * ((-10.0d0 + 9.0d0 * (xi1**2.0d0 + xi2**2.0d0)) + 18.0d0 * xi2 * (1 + xi2))
          dNdxi(11,2) = 9.0d0/32.0d0 * dw4 * (-2.0d0 * xi2 * dw5 + 3.0d0 * (1- xi2 * xi2))
          dNdxi(12,2) = 9.0d0/32.0d0 * dw4 * (-2.0d0 * xi2 * dw6 + 3.0d0 * (1- xi2 * xi2))

!    illegal interpolation code
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
          dNdxi(1,1) = -0.125d0*(1-xi2)*(1-xi3)
          dNdxi(1,2) = -0.125d0*(1-xi1)*(1-xi3)
          dNdxi(1,3) = -0.125d0*(1-xi1)*(1-xi2)
          dNdxi(2,1) =  0.125d0*(1-xi2)*(1-xi3)
          dNdxi(2,2) = -0.125d0*(1+xi1)*(1-xi3)
          dNdxi(2,3) = -0.125d0*(1+xi1)*(1-xi2)
          dNdxi(3,1) =  0.125d0*(1+xi2)*(1-xi3)
          dNdxi(3,2) =  0.125d0*(1+xi1)*(1-xi3)
          dNdxi(3,3) = -0.125d0*(1+xi1)*(1+xi2)
          dNdxi(4,1) = -0.125d0*(1+xi2)*(1-xi3)
          dNdxi(4,2) =  0.125d0*(1-xi1)*(1-xi3)
          dNdxi(4,3) = -0.125d0*(1-xi1)*(1+xi2)
          dNdxi(5,1) = -0.125d0*(1-xi2)*(1+xi3)
          dNdxi(5,2) = -0.125d0*(1-xi1)*(1+xi3)
          dNdxi(5,3) =  0.125d0*(1-xi1)*(1-xi2)
          dNdxi(6,1) =  0.125d0*(1-xi2)*(1+xi3)
          dNdxi(6,2) = -0.125d0*(1+xi1)*(1+xi3)
          dNdxi(6,3) =  0.125d0*(1+xi1)*(1-xi2)
          dNdxi(7,1) =  0.125d0*(1+xi2)*(1+xi3)
          dNdxi(7,2) =  0.125d0*(1+xi1)*(1+xi3)
          dNdxi(7,3) =  0.125d0*(1+xi1)*(1+xi2)
          dNdxi(8,1) = -0.125d0*(1+xi2)*(1+xi3)
          dNdxi(8,2) =  0.125d0*(1-xi1)*(1+xi3)
          dNdxi(8,3) =  0.125d0*(1-xi1)*(1+xi2)

!       20 nodes
        elseif (iplc .eq. 'v20') then
          dNdxi(1,1) =  0.125d0 *(1-xi2)*(1-xi3)*(1+2*xi1+xi2+xi3)
          dNdxi(1,2) =  0.125d0 *(1-xi1)*(1-xi3)*(1+2*xi2+xi1+xi3)
          dNdxi(1,3) =  0.125d0 *(1-xi1)*(1-xi2)*(1+2*xi3+xi1+xi2)
          dNdxi(2,1) = -0.5d0   *(1-xi2)*(1-xi3)* xi1
          dNdxi(2,2) = -0.25d0  *(1-xi3)*(1-xi1**2)
          dNdxi(2,3) = -0.25d0 *(1-xi2)*(1-xi1**2)
          dNdxi(3,1) =  0.125d0 *(1-xi2)*(1-xi3)*(-1+2*xi1-xi2-xi3)
          dNdxi(3,2) =  0.125d0 *(1+xi1)*(1-xi3)*(1+2*xi2-xi1+xi3)
          dNdxi(3,3) =  0.125d0 *(1+xi1)*(1-xi2)*(1+2*xi3-xi1+xi2)
          dNdxi(4,1) =  0.25d0  *(1-xi3)*(1-xi2**2)
          dNdxi(4,2) = -0.5d0   *(1+xi1)*(1-xi3)* xi2
          dNdxi(4,3) = -0.25d0  *(1+xi1)*(1-xi2**2)
          dNdxi(5,1) =  0.125d0 *(1+xi2)*(1-xi3)*(-1+2*xi1+xi2-xi3)
          dNdxi(5,2) =  0.125d0 *(1+xi1)*(1-xi3)*(-1+2*xi2+xi1-xi3)
          dNdxi(5,3) =  0.125d0 *(1+xi1)*(1+xi2)*(1+2*xi3-xi1-xi2)
          dNdxi(6,1) = -0.5d0   *(1+xi2)*(1-xi3)* xi1
          dNdxi(6,2) =  0.25d0  *(1-xi3)*(1-xi1**2)
          dNdxi(6,3) = -0.25d0  *(1+xi2)*(1-xi1**2)
          dNdxi(7,1) =  0.125d0 *(1+xi2)*(1-xi3)*(1+2*xi1-xi2+xi3)
          dNdxi(7,2) =  0.125d0 *(1-xi1)*(1-xi3)*(-1+2*xi2-xi1-xi3)
          dNdxi(7,3) =  0.125d0 *(1-xi1)*(1+xi2)*(1+2*xi3+xi1-xi2)
          dNdxi(8,1) = -0.25d0  *(1-xi3)*(1-xi2**2)
          dNdxi(8,2) = -0.5d0   *(1-xi1)*(1-xi3)* xi2
          dNdxi(8,3) = -0.25d0  *(1-xi1)*(1-xi2**2)
          dNdxi(9,1) = -0.25d0  *(1-xi2)*(1-xi3**2)
          dNdxi(9,2) = -0.25d0  *(1-xi1)*(1-xi3**2)
          dNdxi(9,3) = -0.5d0   *(1-xi1)*(1-xi2)* xi3
          dNdxi(10,1) =  0.25d0  *(1-xi2)*(1-xi3**2)
          dNdxi(10,2) = -0.25d0  *(1+xi1)*(1-xi3**2)
          dNdxi(10,3) = -0.5d0   *(1+xi1)*(1-xi2)* xi3
          dNdxi(11,1) =  0.25d0  *(1+xi2)*(1-xi3**2)
          dNdxi(11,2) =  0.25d0  *(1+xi1)*(1-xi3**2)
          dNdxi(11,3) = -0.5d0   *(1+xi1)*(1+xi2)* xi3
          dNdxi(12,1) = -0.25d0  *(1+xi2)*(1-xi3**2)
          dNdxi(12,2) =  0.25d0  *(1-xi1)*(1-xi3**2)
          dNdxi(12,3) = -0.5d0   *(1-xi1)*(1+xi2)* xi3
          dNdxi(13,1) =  0.125d0 *(1-xi2)*(1+xi3)*(1+2*xi1+xi2-xi3)
          dNdxi(13,2) =  0.125d0 *(1-xi1)*(1+xi3)*(1+2*xi2+xi1-xi3)
          dNdxi(13,3) =  0.125d0 *(1-xi1)*(1-xi2)*(-1+2*xi3-xi1-xi2)
          dNdxi(14,1) = -0.5d0   *(1-xi2)*(1+xi3)* xi1
          dNdxi(14,2) = -0.25d0  *(1+xi3)*(1-xi1**2)
          dNdxi(14,3) =  0.25d0  *(1-xi2)*(1-xi1**2)
          dNdxi(15,1) =  0.125d0 *(1-xi2)*(1+xi3)*(-1+2*xi1-xi2+xi3)
          dNdxi(15,2) =  0.125d0 *(1+xi1)*(1+xi3)*(1+2*xi2-xi1-xi3)
          dNdxi(15,3) =  0.125d0 *(1+xi1)*(1-xi2)*(-1+2*xi3+xi1-xi2)
          dNdxi(16,1) =  0.25d0   *(1+xi3)*(1-xi2**2)
          dNdxi(16,2) = -0.5d0   *(1+xi1)*(1+xi3)* xi2
          dNdxi(16,3) =  0.25d0  *(1+xi1)*(1-xi2**2)
          dNdxi(17,1) =  0.125d0 *(1+xi2)*(1+xi3)*(-1+2*xi1+xi2+xi3)
          dNdxi(17,2) =  0.125d0 *(1+xi1)*(1+xi3)*(-1+2*xi2+xi1+xi3)
          dNdxi(17,3) =  0.125d0 *(1+xi1)*(1+xi2)*(-1+2*xi3+xi1+xi2)
          dNdxi(18,1) = -0.5d0   *(1+xi2)*(1+xi3)* xi1
          dNdxi(18,2) =  0.25d0  *(1+xi3)*(1-xi1**2)
          dNdxi(18,3) =  0.25d0  *(1+xi2)*(1-xi1**2)
          dNdxi(19,1) =  0.125d0 *(1+xi2)*(1+xi3)*(1+2*xi1-xi2-xi3)
          dNdxi(19,2) =  0.125d0 *(1-xi1)*(1+xi3)*(-1+2*xi2-xi1+xi3)
          dNdxi(19,3) =  0.125d0 *(1-xi1)*(1+xi2)*(-1+2*xi3-xi1+xi2)
          dNdxi(20,1) = -0.25d0  *(1+xi3)*(1-xi2**2)
          dNdxi(20,2) = -0.5d0   *(1-xi1)*(1+xi3)* xi2
          dNdxi(20,3) =  0.25d0  *(1-xi1)*(1-xi2**2)

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
          dNdxi(1,1) = 1.0d0
          dNdxi(1,2) = 0.0d0
          dNdxi(1,3) = 0.0d0
          dNdxi(2,1) = 0.0d0
          dNdxi(2,2) = 1.0d0
          dNdxi(2,3) = 0.0d0
          dNdxi(3,1) = 0.0d0
          dNdxi(3,2) = 0.0d0
          dNdxi(3,3) = 1.0d0
          dNdxi(4,1) = -1.0d0
          dNdxi(4,2) = -1.0d0
          dNdxi(4,3) = -1.0d0

!       10 nodes
        elseif (iplc .eq. 'h10') then
          dNdxi(1,1) = 4.0d0 * xi1 - 1.0d0
          dNdxi(1,2) = 0.0d0
          dNdxi(1,3) = 0.0d0
          dNdxi(2,1) = 0.0d0
          dNdxi(2,2) = 4.0d0 * xi2 - 1.0d0
          dNdxi(2,3) = 0.0d0
          dNdxi(3,1) = 0.0d0
          dNdxi(3,2) = 0.0d0
          dNdxi(3,3) = 4.0d0 * xi3 - 1.0d0
          dNdxi(4,1) = 1.0d0 - 4.0d0 * xi4
          dNdxi(4,2) = 1.0d0 - 4.0d0 * xi4
          dNdxi(4,3) = 1.0d0 - 4.0d0 * xi4
          dNdxi(5,1) = 4.0d0 * xi2
          dNdxi(5,2) = 4.0d0 * xi1
          dNdxi(5,3) = 0.0d0
          dNdxi(6,1) = 0.0d0
          dNdxi(6,2) = 4.0d0 * xi3
          dNdxi(6,3) = 4.0d0 * xi2
          dNdxi(7,1) = 4.0d0 * xi3
          dNdxi(7,2) = 0.0d0
          dNdxi(7,3) = 4.0d0 * xi1
          dNdxi(8,1) = 4.0d0 * (xi4-xi1)
          dNdxi(8,2) = -4.0d0*xi1
          dNdxi(8,3) = -4.0d0*xi1
          dNdxi(9,1) = -4.0d0*xi2
          dNdxi(9,2) = 4.0d0 * (xi4-xi2)
          dNdxi(9,3) = -4.0d0*xi2
          dNdxi(10,1)= -4.0d0*xi3
          dNdxi(10,2)= -4.0d0*xi3
          dNdxi(10,3)= 4.0d0 * (xi4-xi3)

!       illegal interpolation code
        else
          write(errmsgtxt, 9990) iplc
 9990     format('Illegal interpolation code: ', a3, '.')
!          call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal interpolation code"
          STOP
        endif



!  Illegal interpolation code
!  --------------------------

      else
        write(errmsgtxt, 9995) iplc
 9995   format('Illegal interpolation code: ', a3, '.')
!        call mexerrmsgtxt(errmsgtxt)
        write(*,*) "Illegal interpolation code"
        STOP
      endif



        END SUBROUTINE gdndxi
END MODULE mod_gdndxi
