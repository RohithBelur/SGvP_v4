MODULE mod_find

    contains

    subroutine find(a, x, b)
    
    ! subroutine similar to find function in matlab
    ! a -> input vector (to be allocated and deallocated in main function)
    ! x -> quantity of intrest
    !
    ! b -> output vector of positions (do not alloate in main function nut deallocate)
    
    implicit none
    
    ! INPUT
    integer*8, allocatable, dimension(:,:), intent(in)       :: a
    integer*8, intent(in)                                    :: x
    
    ! OUTPUT
    integer*8, allocatable, dimension(:,:), intent(out)      :: b
    
    ! WORKSPACE
    integer*8, allocatable, dimension(:,:)                   :: loc
    integer*8                                                :: i, j, k
    
    character(len=100)                                       :: txt
    ! -------------------------------------------------------------------
    
!     write(txt,'(i10)') size(a,1)
!     call mexPrintf(txt)
!     call mexPrintf('\n')
    
    if (allocated(loc)) deallocate(loc)
    allocate (loc(size(a,1),1))

    j = 1
    if (size(a,1) .EQ. 1) then
        do i = 1,size(a,2)
            if (a(i,1) .EQ. x) then
                loc(j,1) = i
                j = j+1
            end if
        end do
    else
        do i = 1,size(a,1)
        do k = 1,size(a,2)
            if (a(i,k) .EQ. x) then
                loc(j,1) = i
                j = j+1
            end if
        end do
        end do
    end if
    
    if (allocated(b)) deallocate(b)
    allocate (b(j-1,1))
    b(:,1) = pack(loc, loc/=0)
    
    deallocate(loc)
!     deallocate(a)
!     deallocate(b)

!     call mexPrintf('end find')
!     call mexPrintf('\n')
    
    end subroutine find

END MODULE mod_find
