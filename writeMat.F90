MODULE write_RealMat

    contains
    
    subroutine write_RealMat(A)
    
    implicit none
    
    real*8, dimension(:,:) :: A
    write(*,*)
    
    do i = lbound(a,1), ubound(a,1)
        write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
    end do
    end subroutine write_RealMat
    
END MODULE write_RealMat

! subroutine write_IntMat(A)
!    integer*8, dimension(:,:) :: A
!    write(*,*)
!    
!    do i = lbound(a,1), ubound(a,1)
!       write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
!    end do
! end subroutine write_IntMat
! 
! subroutine write_vector(A)
!    real*8, dimension(:) :: A
!    write(*,*)
!    
!    do i = 1, ubound(a,1)
!       write(*,*) (a(i))
!    end do
! end subroutine write_vector
