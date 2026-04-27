#include "fintrf.h"
!#include "MEXPTR.h"

#ifndef mwsize
#define mwsize integer*4
#endif

#ifndef mwSize
#define mwSize mwsize
#endif

#ifndef mwpointer
#define mwpointer integer*4
#endif

#ifndef mwPointer
#define mwPointer mwpointer
#endif

!#include "mxCopyPtrToCharacter2.for"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Matlab Gate function for Routine Element plnn2t15
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

use mod_plbr2t15Fct

implicit none




! MEX INPUT/OUTPUT --------------------------------------------------------------------------------
    integer*4                                                           :: nlhs, nrhs                 				! ARGUMENT NUMBER
    integer*8                                                           :: plhs(*), prhs(*)                			! ARGUMENT VALUE
    integer*8								:: ip
    integer*8                                                           :: M, N , Mo, No           				! VARIABLE SIZE

! DATA FOR COMPUTATIONAL ROUTINE ------------------------------------------------------------------

  ! INPUT -------------------------------------------------------------------------
    character(len=10)                                                   :: key
    real*8, 	dimension(6,3)                                          :: nodee
    real*8                                                              :: geome
    real*8, 	dimension(1,16)                                         :: matre
    integer*8, dimension(3)                                             :: datae
    real*8, dimension(48*4)                                             :: eldt
    real*8, 	dimension(15,1)                                         :: ue, due
    real*8                                                              :: dt
    real*8                                                              :: incriT, DincriT

! OUTPUT ------------------------------------------------------------------------
    
    real*8,  dimension (:,:), allocatable                               :: out1
    real*8,  dimension (:,:), allocatable                               :: out2
    real*8, dimension (:,:), allocatable                                :: out3

! FUNCTIONS ---------------------------------------------------------------------------------------
    integer*8                                                           :: mxGetPr
    integer*8                                                           :: mxCreateNumericArray
    integer*8                                                           :: mxCreateDoubleMatrix
    integer*8                                                           :: mxIsNumeric
    integer*8                                                           :: mxGetM, mxGetN
    integer*4								:: mxGetClassID
    integer*4                                                           :: mxGetString
    integer*4                                                           :: flag
    
    character(len=10) :: txt1,txt

! CHECK INPUT/OUTPUT
! =================================================================================================
    !if (nrhs.ne.8)  call mexErrMsgTxt('XXX WRONG INPUT NUMBER (8 REQUIRED)\n')                      ! CHECK INPUT COUNT
    !if (nlhs.ne.3)  call mexErrMsgTxt('XXX WRONG OUTPUT NUMBER (3 REQUIRED)\n')                     ! CHECK OUTPUT COUNT

! AGRUMENT HANDLE
! =================================================================================================

! INPUT ------------------------------------------------------------------------------------------------
    

if (nrhs >= 1) then
! Key   
    ip              = mxGetPr(prhs(1))
    M               = mxGetM(prhs(1))
    N               = mxGetN(prhs(1))
    if              (M.NE.1)                        	call mexErrMsgTxt('XXX key is not a line character\n')
    flag = mxGetString(prhs(1), key, 10)
    if (nrhs >=2) then    
    ! Nodes of element
        ip              = mxGetPr(prhs(2))
        M               = mxGetM(prhs(2))
        N               = mxGetN(prhs(2))
        !if 			(( M.ne.6 ).or.( N.ne.3 )) 				call mexErrMsgTxt('XXX nodee IS NOT [ 6 x 3 ]\n')		! CHECK SIZE
        if 			( mxGetClassID(prhs(2)).ne.6 )			call mexErrMsgTxt('XXX nodee IS NOT REAL*8\n')		! CHECK CLASS
        call          mxCopyPtrToReal8( ip ,nodee,M*N)
        if (nrhs >= 3) then    
        ! Geometry of element
            ip              = mxGetPr(prhs(3))
            M               = mxGetM(prhs(3))
            N               = mxGetN(prhs(3))
            if 			  (( M.ne.1 ).or.( N.ne.1 )) 				call mexErrMsgTxt('XXX geome IS NOT [ 1 x 1 ]\n')		! CHECK SIZE
            if 			 ( mxGetClassID(prhs(3)).ne.6 )			call mexErrMsgTxt('XXX geome IS NOT REAL*8\n')				! CHECK CLASS
            call            mxCopyPtrToReal8( ip ,geome,M*N)

            if (nrhs >= 4) then
            ! Property of element
                ip              = mxGetPr(prhs(4))
                M               = mxGetM(prhs(4))
                N               = mxGetN(prhs(4))
                if 			  (( M.ne.1 ).or.( N.ne.16 )) 			call mexErrMsgTxt('XXX matre IS NOT [ 1 x 15 ]\n')			! CHECK SIZE
                call            mxCopyPtrToReal8( ip ,matre,M*N)
    
                if (nrhs >=5) then
                ! Data [Nb Gauss Point, Plane, Punic]
                    ip              = mxGetPr(prhs(5))
                    M               = mxGetM(prhs(5))
                    N               = mxGetN(prhs(5))    
                    if 			  (( M.ne.1 ).or.( N.ne.3 )) 				call mexErrMsgTxt('XXX datae IS NOT [ 1 x 3]\n')			! CHECK SIZE
                    call            mxCopyPtrToInteger8( ip ,datae,M*N)
   
                    if (nrhs >= 6) then
                    ! ELDTE
                        ip              = mxGetPr(prhs(6))
                        M               = mxGetM(prhs(6))
                        N               = mxGetN(prhs(6))
!                         allocate(eldt(N))
                        call            mxCopyPtrToReal8( ip ,eldt,M*N)

                        if (nrhs >= 7) then
                        ! ue
                            ip              = mxGetPr(prhs(7))
                            M               = mxGetM(prhs(7))
                            N               = mxGetN(prhs(7))    
                            !if 			(( M.ne.15 ).or.( N.ne.1 )) 			call mexErrMsgTxt('XXX ue IS NOT [ 15 x 1 ]\n')		! CHECK SIZE
                            if 			( mxGetClassID(prhs(7)).ne.6 )			call mexErrMsgTxt('XXX ue IS NOT REAL*8\n')			! CHECK CLASS
                            call            mxCopyPtrToReal8( ip ,ue,M*N)

                            if (nrhs >= 8) then
                                ! due
                                ip              = mxGetPr(prhs(8))
                                M               = mxGetM(prhs(8))
                                N               = mxGetN(prhs(8))    
                                !if 			(( M.ne.15 ).or.( N.ne.1 )) 			call mexErrMsgTxt('XXX ve IS NOT [ 15 x 1 ]\n')		! CHECK SIZE
                                if 			( mxGetClassID(prhs(8)).ne.6 )			call mexErrMsgTxt('XXX due IS NOT REAL*8\n')			! CHECK CLASS
                                call            mxCopyPtrToReal8( ip ,due,M*N)
                                
                                if (nrhs >= 9) then
                                    ! dt
                                    ip = mxGetPr(prhs(9))
                                    M  = mxGetM(prhs(9))
                                    N  = mxGetN(prhs(9))
                                    if (mxGetClassID(prhs(9)).ne.6) call mexErrMsgTxt('XXX dt IS NOT REAL*8\n')
                                    call mxCopyPtrToReal8(ip, dt, M*N)
                                    
                                    if (nrhs >= 10) then
                                    ! incriT
                                    ip              = mxGetPr(prhs(10))
                                    M               = mxGetM(prhs(10))
                                    N               = mxGetN(prhs(10))
                                    if ( (M .ne. 1) .or. (N .ne. 1))    call mexErrMsgTxt('XXX incriT IS NOT [ 1 x 1 ]\n')		! CHECK SIZE
                                    if ( mxGetClassID(prhs(10)) .ne. 6)  call mexErrMsgTxt('XXX incriT IS NOT REAL*8\n')	        ! CHECK CLASS
                                    call            mxCopyPtrToReal8(ip, incriT, M*N)
                                    
                                        if (nrhs >= 11) then
                                        ! DincriT
                                        ip              = mxGetPr(prhs(11))
                                        M               = mxGetM(prhs(11))
                                        N               = mxGetN(prhs(11))
                                        if ( (M .ne. 1) .or. (N .ne. 1))    call mexErrMsgTxt('XXX DincriT IS NOT [ 1 x 1 ]\n')		! CHECK SIZE
                                        if ( mxGetClassID(prhs(11)) .ne. 6) call mexErrMsgTxt('XXX DincriT IS NOT REAL*8\n')	        ! CHECK CLASS
                                        call            mxCopyPtrToReal8(ip, DincriT, M*N)
                                    
                                        else
                                        DincriT = 0.d0
                                        end if
                                    else
                                    DincriT = 0.d0
                                    incriT = 0.d0
                                    end if
                                else
                                    DincriT = 0.d0
                                    incriT = 0.d0
                                    due=0.0d0
                                endif
                            else
                                DincriT = 0.d0
                                incriT = 0.d0
                                ue=0.0d0
                                due=0.0d0
                            endif
                        else
                            DincriT = 0.d0
                            incriT = 0.d0
    !                         allocate(eldt(datae(1)*42))
                            ue=0.0d0
                            due=0.0d0
                            eldt=0.0d0
                        endif
                    else
                        DincriT = 0.d0
                        incriT = 0.d0
                        datae=0
    !                     allocate(eldt(datae(1)*42))
                        ue=0.0d0
                        due=0.0d0
                        eldt=0.0d0
                    endif
                else
                    DincriT = 0.d0
                    incriT = 0.d0
                    matre=0.0
                    datae=0
    !                 allocate(eldt(datae(1)*42))
                    ue=0.0d0
                    due=0.0d0
                    eldt=0.0d0
                endif
            else
                DincriT = 0.d0
                incriT = 0.d0
                matre=0.0
                datae=0
!                 allocate(eldt(datae(1)*42))
                ue=0.0d0
                due=0.0d0
                eldt=0.0d0
                dt = 0.0d0
            endif
        else
            DincriT = 0.d0
            incriT = 0.d0
            geome=0.0d0
            matre=0.0
            datae=0
!             allocate(eldt(datae(1)*42))
            ue=0.0d0
            due=0.0d0  
            eldt=0.0d0
            dt = 0.0d0
        endif
    else
        DincriT = 0.d0
        incriT = 0.d0
        nodee=0.0d0
        geome=0.0d0
        matre=0.0
        datae=0
!         allocate(eldt(datae(1)*42))
        ue=0.0d0
        due=0.0d0 
        eldt=0.0d0
        dt = 0.0d0
    endif
else
    call mexErrMsgTxt('Illegal number of inputs \n')
endif

! CALL COMPUTATIONAL ROUTINE
! =================================================================================================
!     call mexPrintf('enter subroutine')
   call plbr2t15Fct(key, nodee, geome, matre, datae, eldt, ue, due, dt, incriT, DincriT, out1, out2, out3)


! FINISH
! =================================================================================================
 
! deallocate(eldt)

if (nlhs >= 1)then	
    Mo=size(out1,1)
   	No=size(out1,2)
    plhs(1)         = mxCreateDoubleMatrix(Mo,No,0)                                              ! INITIALIZE MATLAB ARRAY    
    ip              = mxGetPr(plhs(1)) 
    call            mxCopyReal8ToPtr(out1,ip,Mo*No)                                         ! EXPORT
    if (nlhs >= 2)then
        Mo=size(out2,1)
   	    No=size(out2,2)
        plhs(2)         = mxCreateDoubleMatrix(Mo,No,0)                                              ! INITIALIZE MATLAB ARRAY
        ip              = mxGetPr(plhs(2))
        call            mxCopyReal8ToPtr(out2,ip, Mo*No)                                         ! EXPORT
        if (nlhs >= 3)then
        Mo=size(out3,1)
           No=size(out3,2)
           plhs(3)         = mxCreateDoubleMatrix(Mo,No,0)                                            ! INITIALIZE MATLAB ARRAY
           ip              = mxGetPr(plhs(3)) 
           call            mxCopyReal8ToPtr(out3,ip,Mo*No)                                         ! EXPORT
        else
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3=0.0d0
        endif
    else
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))

        out2=0.0d0
        out3=0.0d0
    endif
else
    call mexErrMsgTxt('Illegal number of outputs \n')  
endif


end subroutine mexFunction
