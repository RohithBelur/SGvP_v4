MODULE mod_gintp

implicit none

CONTAINS

      SUBROUTINE gintp(iplc, nintp, xi, w)
!
!.....FORTRAN version of gintp.m
!...
!...  subroutine gintp(iplc, nintp, xi, ldxi, w)
!...
!...  Generation of integration point information: local coordinates and
!...  weight factors for numerical integration.
!...
!...  input:
!...    iplc     :  interpolation code
!...    nintp    :  number of integration points
!...    ldxi     :  leading dimension of xi-array
!...
!...  output:
!...    xi       :  local coordinates of the integration points
!...    w        :  corresponding weight factors
!...
!
!.....Variable declaration

!
! INPUT -------------------------------------------------------------------------
    character(len=3),intent(in)                         :: iplc
    integer*8,intent(in)                                :: nintp

! OUTPUT ------------------------------------------------------------------------
    real*8,intent(out)                                  :: xi(nintp, *), w(*)
!
! Workspace ---------------------------------------------------------------------
    character(len=80)                                   :: errmsgtxt
    real*8                                              :: dw1, dw2, dw3, dw4

! External routines -------------------------------------------------------------

    !MATLAB matrix manipulation
   ! external          mexerrmsgtxt


! Executable statements ---------------------------------------------------------

! One-dimensional
! ---------------
!

      if (iplc(:1) .eq. 'b') then

!       1 integration point
        if (nintp .eq. 1) then
          xi(1, 1) =  0.0d0
          w(1)     =  2.0d0

!       2 integraton points
        elseif (nintp .eq. 2) then
          dw1      =  1.0d0 / sqrt(3.0d0)
          xi(1, 1) = -dw1
          xi(2, 1) =  dw1
          w(1)     =  1.0d0
          w(2)     =  1.0d0

!       3 integraton points
        elseif (nintp .eq. 3) then
          dw1      =  sqrt(3.0d0/5.0d0)
          xi(1, 1) = -dw1
          xi(2, 1) =  0.0d0
          xi(3, 1) =  dw1
          dw1      =  5.0d0/9.0d0
          dw2      =  8.0d0/9.0d0
          w(1)     =  dw1
          w(2)     =  dw2
          w(3)     =  dw1

!       4 integraton points
        elseif (nintp .eq. 4) then
          xi(1,1)   = -0.861136311594953d0
          xi(2,1)   = -0.339981043584856d0
          xi(3,1)   =  0.339981043584856d0
          xi(4,1)   =  0.861136311594953d0
          w(1)      =  0.347854845137454d0
          w(2)      =  0.652145154862546d0
          w(3)      =  0.652145154862546d0
          w(4)      =  0.347854845137454d0

!       illegal number of integration points
        else
          write(errmsgtxt, 9960) nintp
 9960     format('Illegal number of integration points: ', i3, '.')
 !         call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal number of integration points"
          STOP
        endif


! Two-dimensional, triangular
! ---------------------------

      elseif (iplc(:1) .eq. 't') then
!       1 integration point
        if (nintp .eq. 1) then
          dw1      =  1.0d0 / 3.0d0
          xi(1, 1) =  dw1
          xi(1, 2) =  dw1
          xi(1, 3) =  dw1
          w(1)     =  1.0d0

!       3 integraton points
        elseif (nintp .eq. 3) then
          xi(1, 1) =  0.5d0
          xi(1, 2) =  0.5d0
          xi(1, 3) =  0.0d0
          xi(2, 1) =  0.0d0
          xi(2, 2) =  0.5d0
          xi(2, 3) =  0.5d0
          xi(3, 1) =  0.5d0
          xi(3, 2) =  0.0d0
          xi(3, 3) =  0.5d0
          dw1      =  1.0d0 / 3.0d0
          w(1)     =  dw1
          w(2)     =  dw1
          w(3)     =  dw1

!       4 integraton points
        elseif (nintp .eq. 4) then
          dw1      =  1.0d0 / 3.0d0
          xi(1, 1:3) =  [dw1, dw1, dw1]
          xi(2, 1:3) =  [0.6d0, 0.2d0, 0.2d0]
          xi(3, 1:3) =  [0.2d0, 0.6d0, 0.2d0]
          xi(4, 1:3) =  [0.2d0, 0.2d0, 0.6d0]
          dw1      = -9.0d0 / 16.0d0
          dw2      =  25.0d0 / 48.0d0
          w(1)     =  dw1
          w(2)     =  dw2
          w(3)     =  dw2
          w(4)     =  dw2

!       6 integraton points   MODIFIER
        elseif (nintp .eq. 6) then
          

          dw1 = 0.816847572980459d0
          dw2 = 0.091576213509771d0
          dw3 = 0.108103018168070d0
          dw4 = 0.445948490915965d0
          xi(1, 1) =  dw1
          xi(1, 2) =  dw2
          xi(1, 3) =  dw2
          xi(2, 1) =  dw2
          xi(2, 2) =  dw1
          xi(2, 3) =  dw2
          xi(3, 1) =  dw2
          xi(3, 2) =  dw2
          xi(3, 3) =  dw1
          xi(4, 1) =  dw3
          xi(4, 2) =  dw4
          xi(4, 3) =  dw4
          xi(5, 1) =  dw4
          xi(5, 2) =  dw3
          xi(5, 3) =  dw4
          xi(6, 1) =  dw4
          xi(6, 2) =  dw4
          xi(6, 3) =  dw3
          dw1      =  0.109951743655322d0
          dw2      =  0.223381589678011d0
          w(1)     =  dw1
          w(2)     =  dw1
          w(3)     =  dw1
          w(4)     =  dw2
          w(5)     =  dw2
          w(6)     =  dw2



!       7 integraton points
        elseif (nintp .eq. 7) then
          dw1      =  1.0d0 / 3.0d0
          xi(1, 1) =  dw1
          xi(1, 2) =  dw1
          xi(1, 3) =  dw1
          dw1      =  0.059715871789770d0
          dw2      =  0.470142064105115d0
          xi(2, 1) =  dw1
          xi(2, 2) =  dw2
          xi(2, 3) =  dw2
          xi(3, 1) =  dw2
          xi(3, 2) =  dw1
          xi(3, 3) =  dw2
          xi(4, 1) =  dw2
          xi(4, 2) =  dw2
          xi(4, 3) =  dw1
          dw1      =  0.797426985353087d0
          dw2      =  0.101286507323456d0
          xi(5, 1) =  dw1
          xi(5, 2) =  dw2
          xi(5, 3) =  dw2
          xi(6, 1) =  dw2
          xi(6, 2) =  dw1
          xi(6, 3) =  dw2
          xi(7, 1) =  dw2
          xi(7, 2) =  dw2
          xi(7, 3) =  dw1
          dw1      =  0.225030000300000d0
          w(1)     =  dw1
          dw1      =  0.132394152788506d0
          w(2)     =  dw1
          w(3)     =  dw1
          w(4)     =  dw1
          dw1      =  0.125939180544827d0
          w(5)     =  dw1
          w(6)     =  dw1
          w(7)     =  dw1

!       illegal number of integration points
        else
          write(errmsgtxt, 9970) nintp
 9970     format('Illegal number of integration points: ', i3, '.')
 !         call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal number of integration points"
          STOP
        endif


! Two-dimensional, quadrilateral
! ------------------------------

      elseif (iplc(:1) .eq. 'q') then

!       1 integration point
        if (nintp .eq. 1) then
          xi(1, 1) =  0.0d0
          xi(1, 2) =  0.0d0
          w(1)     =  4.0d0

!       4 integration points
        elseif (nintp .eq. 4) then
          dw1      =  1.0d0 / sqrt(3.0d0)
          xi(1, 1) = -dw1
          xi(1, 2) = -dw1
          xi(2, 1) =  dw1
          xi(2, 2) = -dw1
          xi(3, 1) =  dw1
          xi(3, 2) =  dw1
          xi(4, 1) = -dw1
          xi(4, 2) =  dw1
          w(1)     =  1.0d0
          w(2)     =  1.0d0
          w(3)     =  1.0d0
          w(4)     =  1.0d0

!       9 integration points
        elseif (nintp .eq. 9) then
          dw1      =  sqrt(3.0d0/5.0d0)
          xi(1, 1) = -dw1
          xi(1, 2) = -dw1
          xi(2, 1) =  0.0d0
          xi(2, 2) = -dw1
          xi(3, 1) =  dw1
          xi(3, 2) = -dw1
          xi(4, 1) =  dw1
          xi(4, 2) =  0.0d0
          xi(5, 1) =  dw1
          xi(5, 2) =  dw1
          xi(6, 1) =  0.0d0
          xi(6, 2) =  dw1
          xi(7, 1) = -dw1
          xi(7, 2) =  dw1
          xi(8, 1) = -dw1
          xi(8, 2) =  0.0d0
          xi(9, 1) =  0.0d0
          xi(9, 2) =  0.0d0
          dw1      =  25.0d0 / 81.0d0
          dw2      =  40.0d0 / 81.0d0
          w(1)     =  dw1
          w(2)     =  dw2
          w(3)     =  dw1
          w(4)     =  dw2
          w(5)     =  dw1
          w(6)     =  dw2
          w(7)     =  dw1
          w(8)     =  dw2
          w(9)     =  64.0d0 / 81.0d0

!       illegal number of integration points
        else
          write(errmsgtxt, 9980) nintp
 9980     format('Illegal number of integration points: ', i3, '.')
 !         call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal number of integration points"
          STOP
        endif


!  Three-dimensional, quadrilateral
!  --------------------------------

      elseif (iplc(:1) .eq. 'v') then

!       1 integration point
        if (nintp .eq. 1) then
          xi(1, 1) =  0.0d0
          xi(1, 2) =  0.0d0
          xi(1, 3) =  0.0d0
          w(1)     =  4.0d0

!       8 integration points
        elseif (nintp .eq. 8) then
          dw1      =  1.0d0 / sqrt(3.0d0)
          xi(1, 1) =  dw1
          xi(1, 2) =  dw1
          xi(1, 3) =  dw1
          xi(2, 1) = -dw1
          xi(2, 2) =  dw1
          xi(2, 3) =  dw1
          xi(3, 1) = -dw1
          xi(3, 2) = -dw1
          xi(3, 3) =  dw1
          xi(4, 1) =  dw1
          xi(4, 2) = -dw1
          xi(4, 3) =  dw1
          xi(5, 1) =  dw1
          xi(5, 2) =  dw1
          xi(5, 3) = -dw1
          xi(6, 1) = -dw1
          xi(6, 2) =  dw1
          xi(6, 3) = -dw1
          xi(7, 1) = -dw1
          xi(7, 2) = -dw1
          xi(7, 3) = -dw1
          xi(8, 1) =  dw1
          xi(8, 2) = -dw1
          xi(8, 3) = -dw1

          w(1)     =  1.0d0
          w(2)     =  1.0d0
          w(3)     =  1.0d0
          w(4)     =  1.0d0
          w(5)     =  1.0d0
          w(6)     =  1.0d0
          w(7)     =  1.0d0
          w(8)     =  1.0d0

!       27 integration points
        elseif (nintp .eq. 27) then
          dw1      =  sqrt(3.0d0/5.0d0)
          xi(1, 1) = -dw1
          xi(1, 2) = -dw1
          xi(1, 3) = -dw1
          xi(2, 1) =  0.0d0
          xi(2, 2) = -dw1
          xi(2, 3) = -dw1
          xi(3, 1) =  dw1
          xi(3, 2) = -dw1
          xi(3, 3) = -dw1
          xi(4, 1) =  dw1
          xi(4, 2) =  0.0d0
          xi(4, 3) = -dw1
          xi(5, 1) =  dw1
          xi(5, 2) =  dw1
          xi(5, 3) = -dw1
          xi(6, 1) =  0.0d0
          xi(6, 2) =  dw1
          xi(6, 3) = -dw1
          xi(7, 1) = -dw1
          xi(7, 2) =  dw1
          xi(7, 3) = -dw1
          xi(8, 1) = -dw1
          xi(8, 2) =  0.0d0
          xi(8, 3) = -dw1
          xi(9, 1) =  0.0d0
          xi(9, 2) =  0.0d0
          xi(9, 3) = -dw1
          xi(10,1) = -dw1
          xi(10,2) = -dw1
          xi(10,3) =  0.0d0
          xi(11,1) =  0.0d0
          xi(11,2) = -dw1
          xi(11,3) =  0.0d0
          xi(12,1) =  dw1
          xi(12,2) = -dw1
          xi(12,3) =  0.0d0
          xi(13,1) =  dw1
          xi(13,2) =  0.0d0
          xi(13,3) =  0.0d0
          xi(14,1) =  dw1
          xi(14,2) =  dw1
          xi(14,3) =  0.0d0
          xi(15,1) =  0.0d0
          xi(15,2) =  dw1
          xi(15,3) =  0.0d0
          xi(16,1) = -dw1
          xi(16,2) =  dw1
          xi(16,3) =  0.0d0
          xi(17,1) = -dw1
          xi(17,2) =  0.0d0
          xi(17,3) =  0.0d0
          xi(18,1) =  0.0d0
          xi(18,2) =  0.0d0
          xi(18,3) =  0.0d0
          xi(19,1) = -dw1
          xi(19,2) = -dw1
          xi(19,3) =  dw1
          xi(20,1) =  0.0d0
          xi(20,2) = -dw1
          xi(20,3) =  dw1
          xi(21,1) =  dw1
          xi(21,2) = -dw1
          xi(21,3) =  dw1
          xi(22,1) =  dw1
          xi(22,2) =  0.0d0
          xi(22,3) =  dw1
          xi(23,1) =  dw1
          xi(23,2) =  dw1
          xi(23,3) =  dw1
          xi(24,1) =  0.0d0
          xi(24,2) =  dw1
          xi(24,3) =  dw1
          xi(25,1) = -dw1
          xi(25,2) =  dw1
          xi(25,3) =  dw1
          xi(26,1) = -dw1
          xi(26,2) =  0.0d0
          xi(26,3) =  dw1
          xi(27,1) =  0.0d0
          xi(27,2) =  0.0d0
          xi(27,3) =  dw1
          dw1      =  125.0d0/729.0d0
          dw2      =  200.0d0/729.0d0
          dw3      =  320.0d0/729.0d0
          dw4      =  512.0d0/729.0d0
          w(1)     =  dw1
          w(2)     =  dw2
          w(3)     =  dw1
          w(4)     =  dw2
          w(5)     =  dw1
          w(6)     =  dw2
          w(7)     =  dw1
          w(8)     =  dw2
          w(9)     =  dw3
          w(10)    =  dw2
          w(11)    =  dw3
          w(12)    =  dw2
          w(13)    =  dw3
          w(14)    =  dw2
          w(15)    =  dw3
          w(16)    =  dw2
          w(17)    =  dw3
          w(18)    =  dw4
          w(19)    =  dw1
          w(20)    =  dw2
          w(21)    =  dw1
          w(22)    =  dw2
          w(23)    =  dw1
          w(24)    =  dw2
          w(25)    =  dw1
          w(26)    =  dw2
          w(27)    =  dw3

!       illegal number of integration points
        else
          write(errmsgtxt, 9985) nintp
 9985     format('Illegal number of integration points: ', i3, '.')
 !         call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal number of integration points"
          STOP
        endif


!  Three-dimensional, tetrahedral
!  --------------------------------

      elseif (iplc(:1) .eq. 'h') then

!       1 integration point
        if (nintp .eq. 1) then
          dw1=1.0d0/4.0d0
          xi(1, 1) =  dw1
          xi(1, 2) =  dw1
          xi(1, 3) =  dw1
          xi(1, 4) =  dw1
          w(1)     =  1.0d0

!       4 integration points
        elseif (nintp .eq. 4) then
          dw1=0.585410196624969d0
          dw2=0.138196601125011d0
          xi(1, 1) =  dw1
          xi(1, 2) =  dw2
          xi(1, 3) =  dw2
          xi(1, 4) =  dw2
          xi(2, 1) =  dw2
          xi(2, 2) =  dw1
          xi(2, 3) =  dw2
          xi(2, 4) =  dw2
          xi(3, 1) =  dw2
          xi(3, 2) =  dw2
          xi(3, 3) =  dw1
          xi(3, 4) =  dw2
          xi(4, 1) =  dw2
          xi(4, 2) =  dw2
          xi(4, 3) =  dw2
          xi(4, 4) =  dw1

          dw1=1.0d0/24.0d0
          w(1)     =  dw1
          w(2)     =  dw1
          w(3)     =  dw1
          w(4)     =  dw1
          w(5)     =  dw1
          w(6)     =  dw1
          w(7)     =  dw1
          w(8)     =  dw1

!       illegal number of integration points
        else
          write(errmsgtxt, 9990) nintp
 9990     format('Illegal number of integration points: ', i3, '.')
 !         call mexerrmsgtxt(errmsgtxt)
          write(*,*) "Illegal number of integration points"
          STOP
        endif



! Illegal interpolation code
! --------------------------

      else
        write(errmsgtxt, 9995) iplc
 9995   format('Illegal interpolation code: ', a2, '.')
 !       call mexerrmsgtxt(errmsgtxt)
        write(*,*) "Illegal interpolation code"
        STOP
      endif




    END SUBROUTINE gintp
END MODULE mod_gintp

