MODULE mod_plbr2t15Fct

implicit none

CONTAINS


SUBROUTINE plbr2t15Fct(key, nodee, geome, matre, datae, eldt, ue, due, dt, incriT, DincriT, out1, out2, out3)

use mod_gdndxi
use mod_gintp
use mod_gn
use mod_findDet
use mod_inv
use mod_ghardening

!USE nag_f_ib, Only : f03aaf !, f06eff, f06paf, f06pmf, f06qhf, f06yaf
!
! Modified fortran element for plbr2t15.m
!
!
! Two-dimensional, 6-noded elasto-viscoplastic triangular element with gradient 
! plasticity (Fleck-Hutchinson theory) based on Borg, et el. formulation (2006).  
! Quatradic displacement interpolation and linear plastic  strain.
! Isotropic material properties
! Plane strain formulation
!
!
!
! BORG 2006
!
!
!            !!!!!     FINITE DEFORMATION FORMULATION   !!!!
! 
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!! Element based on the incremental theory !!!!!!!!!!!!!!!!!!
! !!!!!!!!!! MTFEM tricked to work as an incremental procedure !!!!!!!!!!!!
! !!!!!!!!!!!!! with only one iteration per calculation step !!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! input:
!   key      :  input key
!   nodee    :  coordinates of element nodes
!   geome    :  geometry data of element: [A]
!   matre    :  material data of element:
!               in case of power law hardeningmodel:[E nu lint  yield  n
!               blank   m_vp \dot\epsilon_0] % \dot\epsilon_0 is the
!               reference strain rate
!               in case of Voce law hardeningmodel: [E nu lint  yield  theta_0 sigma_v m_vp \dot\epsilon_0]
!   datae    :  element data: [ nintp      hardening model]
!   histe    :  element deformation history at Gauss points:
!   ue       :  relevant components of solution
!   due      :  relevant components of increments
!   dt       :  time increments
!
! input:
!   key      :  input key
!   nodee    :  coordinates of element nodes
!   geome    :  geometry data of element: [A]
!   matre    :  material data of element:
!   datae    :  element data: [nintp]
!   histe    :  element deformation history at Gauss points:
!
!                      <<<<< Previous Parameters >>>>>
!   -------------------------------------------------------------------------------------   
!               [1:3    epsi_prv        total strain (previous step)
!                4      epsipl_prv      effective plastic strain (previous step)
!                5:7    sigma_prv       stress in the directions 11, 22 and 12 (previous step)
!                8      sigma33_prv     out of plane stress component (previous step)
!                9      Ep_prv          accumulated generalised effective 
!                                       plastic strain (previous step)
!                10     Q_prv           generalized effective stress (previous step)
!                11:12  tau_prv         higher order stress (previous step)
!                                       associated to the gradient of the
!                                       effective plastic strain (scalar)
!
!                      <<<<< Current Parameters >>>>>
!   -------------------------------------------------------------------------------------   
!                13:15  epsi            total strain (current step)
!                16     epsipl          effective plastic strain (current step)
!                17:19  sigma           stress in the directions 11, 22 and 12(current step)
!                20     sigma33         out of plane stress component (current step)
!                21     Ep              accumulated generalised effective 
!                                       plastic strain (current step)
!                22     Q               generalized effective stress
!                23:24  tau             higher order stress (current step)
!
!                25     depsipl         effective plastic strain rate FROM Gauss points (previous step)
!                26:27  ddepsipldx      plastic strain spatial gradient rate FROM Gauss points (previous step)
!                28     dg_dEp          hardenning modulous 
!                29:32  m               gradient to the yield surface
!                                       ~df/dsigma
!                33     RESERVED        for "unloading flag" (in developped version)
!                34     ssl_flag        second step loading flag - first
!                                       step set as elastic  to avoid problems with Eqs. 10, 11
!                                       and the remaining ones are viscoplastic
!                35     sigmac          effective stress used to calculate generalised
!                                       effective plastic strain rate
!                                       (viscoplastic model) Eq. 12
!                36     sigmae          vonMises equivalent stress
!                37     dEp             Generalised effective plastic
!                                       strain rate Eq. 13
!                38     gEp             current yield strength Eq. 14             
!                39:40  rho             Kirchhoff higher order stress Eq. 6
!                41:44  veps_pl         plastic strain vector
!                45:48  epsthv          thermal strain vector
!                ]
!   ue       :  relevant components of solution
!   due      :  relevant components of increments
!
! output:
!   out1     :  dependent on 'key':
!               key = 'nnodee'   :  number of element nodes
!                     'ndofne'   :  number of nodal degrees of freedom
!                     'Ke'       :  element stiffness matrix
!                     'finte'    :  stress and plastic strain update
!                      'fvp'       viscoplastic force increment Eq. 27,
!                      28
!                     'histe'    :  processed deformation history
!                     'displa'   :  displacement
!                     'strain'   :  strain
!                     'plshst'   :  plastic strain
!                     'stress'   :  stress
!   out2     :  dependent on 'key':
!               key = 'finte'    :  updated element deformation history
!
!
implicit none
!
! INPUT -------------------------------------------------------------------------
    character(len=*), intent(in)                                        :: key
    real*8, intent(in)                                                  :: nodee(6,3)
    real*8, intent(in)                                                  :: geome
    real*8, intent(in)                                                  :: matre(1,16)
    integer*8, intent(in)                                               :: datae(3)
    real*8, intent(in)                                                  :: eldt(55*datae(1))
    real*8, intent(in)                                                  :: ue(15,1), due(15,1)
    real*8, intent(in)                                                  :: dt
    real*8, intent(in)                                                  :: incriT, DincriT

! OUTPUT ------------------------------------------------------------------------
    real*8, dimension (:,:), allocatable, intent(out)           :: out1
    real*8, dimension (:,:), allocatable, intent(out)           :: out2
    real*8, dimension (:,:), allocatable, intent(out)           :: out3
!
! Workspace ---------------------------------------------------------------------
    real*8, dimension(:), allocatable                           :: Tempeldt
    integer*8                                                   :: i, nintp, iintp, MAXNHISTI, nnode, ii
    character(len=3)                                            :: ipole, ipolu
    integer*8                                                   :: ndof, ndofu, ndofe
    integer*8                                                   :: uix(12), eix(3)
    real*8                                                      :: coeff, tole, k
    real*8, dimension(datae(1),3)                               :: xi
    real*8, dimension(datae(1))                                 :: w
    real*8                                                      :: dNdxi(6,2),dNdx(6,2),dNedxi(3,2),dNedx(3,2)
    real*8                                                      :: displa(2,6), ddispla(2,6)
    real*8, dimension(:,:), allocatable                         :: Iden, I4, I3
    real*8                                                      :: J(2,2), J0(2,2), detJ, invJ(2,2), invJ0(2,2)
    real*8                                                      :: t, E, nu, l_int, Y, theta0, sigma_v, Ei, hi, alphav, gamma, twan
    real*8                                                      :: mvp, deps0, n1
    integer                                                     :: hardening_law
    real*8, dimension(:,:), allocatable                         :: Kuu, Kup, Kpu, Kpp
    real*8, dimension(:,:), allocatable                         :: fu, fe, fte, Et
    real*8                                                      :: H(4,4), Dg(4,4), D11(3,3), D12(3,1), D22, Ht(3,3)
    real*8                                                      :: N(6), Ne(1,3), NeT(3,1), Tr(1,3), trNe(3,1)
    real*8, dimension(:,:), allocatable                         :: ones
    real*8, parameter                                           :: pi = 4.d0 * atan(1.0_8)

    real*8                                                      :: Eu(4,12), Bu(3,12), Be(2,3), Ce(15,15), Ke(15,15)
    real*8                                                      :: m(1,4), h_Ep, pu(1,4), SIGMA1(2,2)
    real*8                                                      :: omega(1,4),OMEGA1(2,2), DS(2,2), ds1(1,4)
    real*8                                                      :: rho(2,1), Q, tau(2,1)
    real*8                                                      :: theta1(2,2), theps(2,2), epsthv(4,1), ag(2,2), dag(2,2)
    real*8                                                      :: delta_epsth(2,2), delta_epsthv(4,1), delta_depsthv(4,1)
    integer                                                     :: ssl_flag, unloading_flag, plastic_flag
    real*8                                                      :: delta_epsi(4,1), epsi(4,1), delta_e(4,1), Ee(2,2)
    real*8                                                      :: depsipl, ddepsipldx(2,1)
    real*8                                                      :: delta_epsipl, epsipl
    real*8                                                      :: delta_depsipl(1,1), delta_ddepsipldx(2,1)
    real*8                                                      :: delta_sigma_j(4,1), delta_sigma(4,1), sigmam, sigma(4,1)
    real*8                                                      :: dev1(4,1), dev(4,1), devs(4,1), sigmae
    real*8                                                      :: Ep, dEp, dg_dEp, gEp, delta_dEp(1,1), delta_Ep, dEp_mem
    real*8                                                      :: sigmac, delta_Q(1,1), delta_rho_c(2,1), delta_rho(2,1)
    real*8                                                      :: delta_tau(2,1), sigmac_prv
    real*8                                                      :: finte(15,1)
    real*8                                                      :: epsi_prv(4,1), epsipl_prv, sigma_prv(4,1), Ep_prv, Q_prv
    real*8                                                      :: tau_prv(2,1)
    real*8                                                      :: gEp_prv, rho_prv(2,1), veps_pl(4,1)
    real*8                                                      :: ad1, delta_depsipl_ad(1,1), xy(1,1)
    
    real*8                                                      :: fvpu(12,1), fvpe(3,1), fvp(15,1), qs(1,1)
    
    real*8, dimension(:,:), allocatable                         :: outn, outi, Ee1, Fe1
    real*8                                                      :: volume, vastrn(1,3), vastrs(1,4), xx(3,3), sigV
    real*8                                                      :: vastpl
    
    real*8                                                      :: vastress(1,4), vastrain(1,3)
    
    character(len=300)                                          :: txt,txt1
    
    ! Logical conditions
    logical :: cond1, cond2, cond3, cond4

    ! Integer flags (0 = false, 1 = true)
    integer :: flag_cond1, flag_cond2, flag_cond3, flag_cond4
    real*8 :: tol_rel_geom, tol_abs_geom, tol_rel_strain, tol_abs_strain
    real*8 :: tol1, tol4, ref1, ref4


! External routines -------------------------------------------------------------

    !MATLAB matrix manipulation
   ! external          mexerrmsgtxt
!NAG
   !   external          f03aaf, f04aaf, f06eff, f06paf, f06pmf, f06qhf, f06yaf

! Executable statements ---------------------------------------------------------

! maximum number of integration point data



!  maximum number of integration point history data
MAXNHISTI = 55
!  number of nodes
nnode = 6
!  indices of the degrees-of-freedom for the displacements
uix = (/1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15/)
!  indices of the degrees-of-freedom for the plastic shear strain
eix = (/3, 8, 13/)
!  total number of degrees-of-freedom for the element
ndof  = 15
!  number of degrees-of-freedom for the displacements
ndofu = 12
!  number of degrees-of-freedom for the plastic strain
ndofe = 3
!  interpolation code for the displacements (quadratic)
ipolu = 't06'
!  interpolation code for the plastic shearing strain (linear)
ipole = 't03'
!  coefficient for constraint on the plastic flow
 coeff = 1e-8
! tolerance for Eq. 13
tole = 1e-20
!
!
! Identity matrix
allocate(Iden(ndofe,ndofe))
Iden(1:ndofe,1:ndofe) = 0.d0
forall(i = 1:ndofe) Iden(i,i) = 1.d0

allocate(I4(3,4))
I4(1:3,1:4) = 0.d0
I4(1,:) = [1.d0,0.d0,0.d0,0.d0]
I4(2,:) = [0.d0,1.d0,0.d0,0.d0]
I4(3,:) = [0.d0,0.d0,1.d0,0.d0]

allocate(I3(4,3))
I3(1:4,1:3) = 0.d0
I3(1,:) = [1.d0,0.d0,0.d0]
I3(2,:) = [0.d0,1.d0,0.d0]
I3(3,:) = [0.d0,0.d0,1.d0]
I3(4,:) = [0.d0,0.d0,0.d0]

! Number of nodes in element
! --------------------------

    if (key.EQ.'nnodee') then
            
      if (allocated(out1)) deallocate(out1)
      allocate(out1(1,1))
      out1 = nnode

      if (allocated(out2)) deallocate(out2)
      allocate(out2(1,1))
      out2=0

      if (allocated(out3)) deallocate(out3)
      allocate(out3(1,1))
      out3=0

! Number of nodal degrees of freedom
! ----------------------------------

    elseif (key.EQ.'ndofne') then

        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,nnode))
        out1(1,1:nnode) = (/3, 2, 3, 2, 3, 2/)

        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2=0

        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3=0
!
!
    elseif (key .eq. 'histe') then
            ! number of integration points
        
            nintp = datae(1)
            
            allocate(Tempeldt(MAXNHISTI*nintp))
            Tempeldt = eldt
!     
            ! initialization
            if (sum(Tempeldt) .eq. 0.d0) then
                
                ! check number of integration points
                if ((nintp .ne. 1) .and. (nintp .ne. 4) .and. (nintp .ne. 9)) then
                    print*,'Illegal number of integration points: ', nintp,'.'
                    STOP
                endif
!             
                ! Initialize history
                Tempeldt(1:(MAXNHISTI*nintp)) = 0.d0
                Ei = matre(1,1)
                hi = Ei
                
                do ii = 1,nintp
                    Tempeldt(MAXNHISTI*(ii-1) + 38) = hi
                end do
            
            ! processing at end of increment
            else
                
                ! copy current values of history parameters to converged ones
                do ii = 1,nintp
                    Tempeldt(MAXNHISTI*(ii-1) + 1)    = Tempeldt(MAXNHISTI*(ii-1) + 13)           ! epsi1
                    Tempeldt(MAXNHISTI*(ii-1) + 2)    = Tempeldt(MAXNHISTI*(ii-1) + 14)           ! epsi2    
                    Tempeldt(MAXNHISTI*(ii-1) + 3)    = Tempeldt(MAXNHISTI*(ii-1) + 15)           ! epsi3
                    Tempeldt(MAXNHISTI*(ii-1) + 4)    = Tempeldt(MAXNHISTI*(ii-1) + 16)           ! epsipl
                    Tempeldt(MAXNHISTI*(ii-1) + 5)    = Tempeldt(MAXNHISTI*(ii-1) + 17)           ! sigma
                    Tempeldt(MAXNHISTI*(ii-1) + 6)    = Tempeldt(MAXNHISTI*(ii-1) + 18)           ! sigma
                    Tempeldt(MAXNHISTI*(ii-1) + 7)    = Tempeldt(MAXNHISTI*(ii-1) + 19)           ! sigma
                    Tempeldt(MAXNHISTI*(ii-1) + 8)    = Tempeldt(MAXNHISTI*(ii-1) + 20)           ! sigma
                    Tempeldt(MAXNHISTI*(ii-1) + 9)    = Tempeldt(MAXNHISTI*(ii-1) + 21)           ! Ep
                    Tempeldt(MAXNHISTI*(ii-1) + 10)   = Tempeldt(MAXNHISTI*(ii-1) + 22)           ! Generalised effective stress Q
                    Tempeldt(MAXNHISTI*(ii-1) + 11)   = Tempeldt(MAXNHISTI*(ii-1) + 23)           ! Higher order stress 1
                    Tempeldt(MAXNHISTI*(ii-1) + 12)   = Tempeldt(MAXNHISTI*(ii-1) + 24)           ! Higher order stress 2
                end do
            
            endif
            
            !   output
            if (allocated(out1)) deallocate(out1)
            allocate(out1(1,(MAXNHISTI*nintp)))
            out1(1,:) = Tempeldt(1:MAXNHISTI*nintp)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
            
            deallocate(Tempeldt)
!
!         
! Element stiffness matrix
! --------------------------------------------
    elseif (key.EQ.'Ke') then
    
        ! thickness
        t = geome
        ! Youngs Modulus
        E = matre(1,1)
        ! Elastic shear modulus
        nu = matre(1,2)
        
        if (incriT .eq. 0.d0) then
        ! rate sensitivity parameter
            mvp = matre(1,4)
        ! Intrensic length
            l_int = matre(1,3)
        else
            mvp = matre(1,15)
        ! Intrensic length
            l_int = matre(1,16)
        end if
        
        ! number of integration points
        nintp = datae(1)

        allocate(Tempeldt(MAXNHISTI*nintp))
        Tempeldt = eldt
        
        ! local coordinates and weight factors of intergration points
        call gintp(ipolu, nintp, xi, w)
    
        ! compute element stiffness matrix in the integration points
        allocate(Kuu(ndofu,ndofu))
        allocate(Kup(ndofu,ndofe))
        allocate(Kpu(ndofe,ndofu))
        allocate(Kpp(ndofe,ndofe))
        
        Kuu(1:ndofu,1:ndofu) = 0.d0
        Kup(1:ndofu,1:ndofe) = 0.d0
        Kpu(1:ndofe,1:ndofu) = 0.d0
        Kpp(1:ndofe,1:ndofe) = 0.d0
        
        ! Hookean matrix plane strain
        H(1:4,1:4) = 0.d0
        H(1,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[1.0d0-nu,  nu,        0.0d0,                  nu      ]
        H(2,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        1.0d0-nu,  0.0d0,                  nu      ]
        H(3,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[0.0d0,     0.0d0,     (1.0d0-2.0d0*nu)/2.0d0, 0.0d0   ]
        H(4,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        nu,        0.0d0,                  1.0d0-nu]
        
        ! initializations
        Eu(1:4,1:ndofu) = 0.d0
        Bu(1:3,1:ndofu) = 0.d0
        Dg(1:4,1:4) = 0.d0
        
        do iintp = 1,nintp
            
            ! read from history
            sigma(1:3,1)       = Tempeldt(MAXNHISTI*(iintp-1)+ 5:7)
            Q                  = Tempeldt(MAXNHISTI*(iintp-1)+ 22)
            rho(1:2,1)         = Tempeldt(MAXNHISTI*(iintp-1)+ 39:40)
            depsipl            = Tempeldt(MAXNHISTI*(iintp-1)+ 25)
            ddepsipldx(1:2,1)  = Tempeldt(MAXNHISTI*(iintp-1)+ 26:27)
            m(1,1:4)           = Tempeldt(MAXNHISTI*(iintp-1)+ 29:32)
            ssl_flag           = Tempeldt(MAXNHISTI*(iintp-1)+ 34)
            sigmac             = Tempeldt(MAXNHISTI*(iintp-1)+ 35)
            dEp                = Tempeldt(MAXNHISTI*(iintp-1)+ 37)
            
            unloading_flag     = Tempeldt(MAXNHISTI*(iintp-1)+ 50)
            
            dEp_mem = Tempeldt(MAXNHISTI*(iintp-1)+ 49)
            
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            ! strain displacement matrix Bu
            dNdx = matmul(dNdxi,invJ)
            
            do ii = 1,6
                Eu(1,2*(ii-1)+1) = dNdx(ii,1)
                Eu(2,2*(ii-1)+2) = dNdx(ii,2)
                Eu(3,2*(ii-1)+1) = dNdx(ii,2)
                Eu(4,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            do ii = 1,6
                Bu(1,2*(ii-1)+1) = dNdx(ii,1)
                Bu(2,2*(ii-1)+2) = dNdx(ii,2)
                Bu(3,2*(ii-1)+1) = dNdx(ii,2)
                Bu(3,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            ! geometric stiffness matrix
            Dg(1,:) = [-sigma(1,1),  0.0d0,        0.0d0,                          -sigma(3,1) ]
            Dg(2,:) = [0.0d0,        -sigma(2,1),  -sigma(3,1),                    0.0d0       ]
            Dg(3,:) = [0.0d0,        -sigma(3,1),  -(sigma(1,1)-sigma(2,1))/2.d0,  -(sigma(2,1)+sigma(1,1))/2.d0 ]
            Dg(4,:) = [-sigma(3,1),  0.0d0,        -(sigma(2,1)+sigma(1,1))/2.d0,   (sigma(1,1)-sigma(2,1))/2.d0 ]
            
            Kuu(1:ndofu,1:ndofu) = Kuu(1:ndofu,1:ndofu) + w(iintp) * t * matmul(matmul(transpose(Bu) , H(1:3,1:3)),Bu)*detJ &
                                                      & + w(iintp) * t * matmul(matmul(transpose(Eu) , Dg) , Eu) * detJ
            
            ! interpoation matrix for plastic strain
            call gdndxi(ipole, xi(iintp, :),nnode/2, dNedxi)
            
            dNedx = matmul(dNedxi,invJ)
            
            Be = transpose(dNedx)
                
            call gn(ipole, xi(iintp, :),Ne)
            NeT = transpose(Ne)
            
            ! account for the contribution of integration points
            if (ssl_flag .eq. 0) then
                ! if elastic or initialization construct Kpp
                k = coeff
                Kpp = Kpp + k*Iden
            else
                ! viscoplastic case
                ! contribution from integration points (eq. 25 and 26 in the paper)
                
                Kpu(1:ndofe,1:ndofu) = Kpu(1:ndofe,1:ndofu) - w(iintp) * t * matmul(matmul(matmul(NeT,m),matmul(H,I3)),Bu)*detJ
                
                if (unloading_flag .eq. 1) then
                    Kpp(1:ndofe,1:ndofe) = 1.d0
                    
                    call mexPrintf('\n')
                    call mexPrintf('Kpp fix !!!')
                    call mexPrintf('\n')
                else
                    Kpp(1:ndofe,1:ndofe) = Kpp(1:ndofe,1:ndofe) &
                            & + w(iintp) * t * (depsipl*(mvp-1.d0)*Q/dEp**2 + sigmac/dEp )*matmul(NeT,Ne) * detJ &
                            & + w(iintp) * t * (l_int/dEp)**2 *(mvp-1.d0)*Q*matmul(NeT,matmul(transpose(ddepsipldx),Be))*detJ &
                            & + w(iintp) * t * (depsipl * (mvp-1.d0)/dEp**2)*matmul(matmul(transpose(Be),rho),Ne)*detJ &
                            & + w(iintp) * t * (l_int/dEp)**2 *(mvp-1.d0)* &
                                                & matmul(matmul(transpose(Be),ddepsipldx),matmul(transpose(rho),Be))*detJ &
                                & + w(iintp) * t * l_int**2 * sigmac/dEp * matmul(transpose(Be),Be)*detJ
                
                endif
            
                

!                 xx = (depsipl*(mvp-1.d0)*Q/dEp**2 + sigmac/dEp )*matmul(NeT,Ne)
                
!                 call mexPrintf('Kpp = ')
!                 write(txt,'(12e15.6)') Kpp(1,1), Kpp(1,2), Kpp(1,3)
!                 call mexPrintf(txt)
!                 call mexPrintf('\n')
!                 
!                 write(txt,'(12e15.6)') Kpp(2,1), Kpp(2,2), Kpp(2,3)
!                 call mexPrintf(txt)
!                 call mexPrintf('\n')
!                 
!                 write(txt,'(12e15.6)') Kpp(3,1), Kpp(3,2), Kpp(3,3)
!                 call mexPrintf(txt)
!                 call mexPrintf('\n')
                
            endif
            
        end do

        ! construct element tangential stiffness matrix
        Ke(1:ndof,1:ndof) = 0.d0
        Ke(uix, uix) = Kuu(1:ndofu,1:ndofu)
        Ke(uix, eix) = Kup(1:ndofu,1:ndofe)
        Ke(eix, uix) = Kpu(1:ndofe,1:ndofu)
        Ke(eix, eix) = Kpp(1:ndofe,1:ndofe)
        
        !   output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(15,15))
        out1 = Ke
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
        deallocate(Kuu)
        deallocate(Kup)
        deallocate(Kpu)
        deallocate(Kpp)
        deallocate(Tempeldt)
        
    ! eigen strain effect on the force vector
    elseif (key .EQ. 'ft') then
        
        allocate(fte(ndof,1))
        fte(1:ndof, 1) = 0.d0
        
        t = geome
        nintp = datae(1)
        E = matre(1,1)
        nu = matre(1,2)
        
        alphav = matre(1,12)
        gamma  = matre(1,13)
        twan   = matre(1,14)
        
        if (incriT .eq. 0.d0) then
            Y = matre(1,6)
        else
            Y = matre(1,6) + (matre(1,7) - matre(1,6))*incriT
        end if
        
        ! local coordinates and weight factors of intergration points
        call gintp(ipolu, nintp, xi, w)

        allocate(Tempeldt(MAXNHISTI*nintp))
        Tempeldt = eldt

        ! initializations
        Eu(1:4,1:ndofu) = 0.d0
        Bu(1:3,1:ndofu) = 0.d0
        H(1:4,1:4) = 0.d0
        
        theta1(1:2,1:2) = 0.d0                              
        ag(1:2,1:2) = 0.d0
        theps(1:2,1:2) = 0.d0
        epsthv(1:4,1) = 0.d0
        dag(1:2,1:2) = 0.d0
        delta_epsth(1:2,1:2) = 0.d0
        delta_epsthv(1:4,1) = 0.d0
        
        do iintp = 1,nintp
            
            ! get current step quantities
            m(1,1:4)           = Tempeldt(MAXNHISTI*(iintp-1)+ 29:32)
            ssl_flag           = Tempeldt(MAXNHISTI*(iintp-1)+ 34)
            
            if (ssl_flag .eq. 0) then
                m(1,1:4) = 0.d0
            end if
            
            ! compute Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            ! strain displacement matrix Bu
            dNdx = matmul(dNdxi,invJ)
            
            do ii = 1,6
                Eu(1,2*(ii-1)+1) = dNdx(ii,1)
                Eu(2,2*(ii-1)+2) = dNdx(ii,2)
                Eu(3,2*(ii-1)+1) = dNdx(ii,2)
                Eu(4,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            do ii = 1,6
                Bu(1,2*(ii-1)+1) = dNdx(ii,1)
                Bu(2,2*(ii-1)+2) = dNdx(ii,2)
                Bu(3,2*(ii-1)+1) = dNdx(ii,2)
                Bu(3,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            
            ! interpoation matrix for plastic strain
            call gdndxi(ipole, xi(iintp, :),nnode/2, dNedxi)
            
            dNedx = matmul(dNedxi,invJ)
            
            Be = transpose(dNedx)
                
            call gn(ipole, xi(iintp, :),Ne)
            NeT = transpose(Ne)
        
            ! Hookean matrix plane strain
            H(1,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[1.0d0-nu,  nu,        0.0d0,                  nu      ]
            H(2,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        1.0d0-nu,  0.0d0,                  nu      ]
            H(3,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[0.0d0,     0.0d0,     (1.0d0-2.0d0*nu)/2.0d0, 0.0d0   ]
            H(4,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        nu,        0.0d0,                  1.0d0-nu]
            
            D11(1:3,1:3) = H(1:3,1:3)
            D12(1:3,1)   = H(1:3,4)
            D22          = H(4,4)
            Ht = D11 - matmul(D12,transpose(D12))/D22
            
            theta1(1,1) =  cos(twan*pi/180.d0)
            theta1(1,2) =  sin(twan*pi/180.d0)
            theta1(2,1) = -sin(twan*pi/180.d0)
            theta1(2,2) =  cos(twan*pi/180.d0)
            
            ag(1,2) = incriT*gamma/2.d0
            ag(2,1) = incriT*gamma/2.d0
            ag(2,2) = incriT*alphav
            
            theps = matmul(matmul(theta1,ag),transpose(theta1))
            
            epsthv(1,1) = theps(1,1)
            epsthv(2,1) = theps(2,2)
            epsthv(3,1) = theps(1,2)*2.d0
            
            dag(1,2) = DincriT*gamma/2.d0
            dag(2,1) = DincriT*gamma/2.d0
            dag(2,2) = DincriT*alphav
            
            delta_epsth = matmul(matmul(theta1,dag),transpose(theta1))
            
            delta_epsthv(1,1) = delta_epsth(1,1)
            delta_epsthv(2,1) = delta_epsth(2,2)
            delta_epsthv(3,1) = delta_epsth(1,2)*2.d0
            
            delta_depsthv = delta_epsthv*dt
            
            fte(uix,1) = fte(uix,1) + w(iintp) * t * matmul(matmul(transpose(Bu),Ht),delta_depsthv(1:3,1))*detJ
            fte(eix,1) = fte(eix,1) - w(iintp) * t * &
                        & matmul(matmul(matmul(transpose(Ne(1:1,1:3)),m(1:1,1:3)),Ht),delta_depsthv(1:3,1))*detJ
            
            
            
        end do
!     
        !   output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(ndof,1))
        out1(1:ndof,1) = fte(1:ndof,1)
        
        deallocate(fte)
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
        deallocate(Tempeldt)
        
        ! stress and plastic strain update
    
    elseif (key.EQ.'finte') then
            
        ! thickness
        t = geome
        ! Youngs Modulus
        E = matre(1,1)
        ! Elastic shear modulus
        nu = matre(1,2)
        
        if (incriT .eq. 0.d0) then
            Y = matre(1,6)
        ! rate sensitivity parameter
            mvp = matre(1,4)
        ! Intrensic length
            l_int = matre(1,3)
        else
            Y = matre(1,6) + (matre(1,7) - matre(1,6))*incriT
        ! rate sensitivity parameter
            mvp = matre(1,15)
        ! Intrensic length
            l_int = matre(1,16)
        end if
        
        ! hardening exponent
        n1 = matre(1,8)
        ! reference strain rate
        deps0 = matre(1,5)
        ! number of integration points
        nintp = datae(1)
        ! hardening law
        hardening_law = datae(2)
        
        alphav = matre(1,12)
        gamma  = matre(1,13)
        twan   = matre(1,14)
        
        ! local coordinates and weight factors of intergration points
        call gintp(ipolu, nintp, xi, w)
        
        ! compute right side in integration points loop
        allocate(fu(ndofu,1))
        fu(1:ndofu,1) = 0.d0    ! displacement dofs
        
        allocate(fe(ndofe,1))
        fe(1:ndofe,1) = 0.d0    ! viscoplastic strain rate dofs

        allocate(Tempeldt(MAXNHISTI*nintp))
        Tempeldt(MAXNHISTI*nintp) = 0.d0
        Tempeldt = eldt
        
        ! initilizations of variables
        SIGMA1(1:2,1:2) = 0.d0
        Eu(1:4,1:ndofu) = 0.d0
        Bu(1:3,1:ndofu) = 0.d0
        delta_epsi(1:4,1) = 0.d0
        epsi(1:4,1) = 0.d0
        OMEGA1(1:2,1:2) = 0.d0
        Ee(1:2,1:2) = 0.d0
        
        delta_depsipl = 0.d0
        delta_ddepsipldx(1:2,1) = 0.d0
        delta_epsipl = 0.d0
        epsipl = 0.d0
        H(1:4,1:4) = 0.d0
        theta1(1:2,1:2) = 0.d0      
        ag(1:2,1:2) = 0.d0
        theps(1:2,1:2) = 0.d0
        epsthv(1:4,1) = 0.d0
        dag(1:2,1:2) = 0.d0
        delta_epsth(1:2,1:2) = 0.d0
        delta_epsthv(1:4,1) = 0.d0
        ds1(1,1:4) = 0.d0
        delta_sigma(1:4,1) = 0.d0
        sigma(1:4,1) = 0.d0
        delta_dEp = 0.d0
        delta_Q = 0.d0
        delta_rho_c(1:2,1) = 0.d0
        delta_rho = 0.d0
        delta_tau = 0.d0
        Q = 0.d0
        rho(1:2,1) = 0.d0
        tau(1:2,1) = 0.d0
        sigmac = 0.d0
        dEp = 0.d0
        delta_Ep = 0.d0
        Ep = 0.d0
        depsipl = 0.d0
        ddepsipldx(1:2,1) = 0.d0
        m(1,1:4) = 0.d0
        
        veps_pl(1:4,1)   = 0.d0
        epsi_prv(1:3,1)  = 0.d0
        epsipl_prv       = 0.d0
        sigma_prv(1:4,1) = 0.d0
        Ep_prv           = 0.d0
        Q_prv            = 0.d0
        tau_prv(1:2,1)   = 0.d0
        gEp_prv          = 0.d0
        rho_prv(1:2,1)   = 0.d0
        
        depsipl           = 0.d0
        ddepsipldx(1:2,1) = 0.d0
        dg_dEp            = 0.d0
        m(1,1:4)          = 0.d0
        ssl_flag          = 0.d0
        sigmac            = 0.d0
        dEp               = 0.d0
        plastic_flag      = 0
        unloading_flag    = 0
        flag_cond1        = 0
        flag_cond2        = 0
        flag_cond3        = 0
        flag_cond4        = 0

        do iintp = 1,nintp
            
            ! read from history the previous increment
            epsi_prv(1:3,1)  = Tempeldt(MAXNHISTI*(iintp-1)+ 1:3)
            epsipl_prv       = Tempeldt(MAXNHISTI*(iintp-1)+ 4)
            sigma_prv(1:4,1) = Tempeldt(MAXNHISTI*(iintp-1)+ 5:8)
            Ep_prv           = Tempeldt(MAXNHISTI*(iintp-1)+ 9)
            Q_prv            = Tempeldt(MAXNHISTI*(iintp-1)+ 10)
            tau_prv(1:2,1)   = Tempeldt(MAXNHISTI*(iintp-1)+ 11:12)
            gEp_prv          = Tempeldt(MAXNHISTI*(iintp-1)+ 38)
            rho_prv(1:2,1)   = Tempeldt(MAXNHISTI*(iintp-1)+ 39:40)
            veps_pl(1:4,1)   = Tempeldt(MAXNHISTI*(iintp-1)+ 41:44)
            
            SIGMA1(1,1) = sigma_prv(1,1)
            SIGMA1(2,2) = sigma_prv(2,1)
            SIGMA1(1,2) = sigma_prv(3,1)
            SIGMA1(2,1) = sigma_prv(3,1)
            
            ! read from history the current values
            depsipl           = Tempeldt(MAXNHISTI*(iintp-1)+ 25)
            ddepsipldx(1:2,1) = Tempeldt(MAXNHISTI*(iintp-1)+ 26:27)
            dg_dEp            = Tempeldt(MAXNHISTI*(iintp-1)+ 28)
            m(1,1:4)          = Tempeldt(MAXNHISTI*(iintp-1)+ 29:32)
            ssl_flag          = Tempeldt(MAXNHISTI*(iintp-1)+ 34)
            sigmac            = Tempeldt(MAXNHISTI*(iintp-1)+ 35)
            dEp               = Tempeldt(MAXNHISTI*(iintp-1)+ 37)
            
            dEp_mem = dEp
            sigmac_prv = sigmac
            
!             call mexPrintf('\n')
!             call mexPrintf('hello')
!             call mexPrintf('\n')

            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
            ddispla = reshape(due(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+(displa-ddispla), dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            ! strain displacement matrix Bu
            dNdx = matmul(dNdxi,invJ)
            
            do ii = 1,6
                Eu(1,2*(ii-1)+1) = dNdx(ii,1)
                Eu(2,2*(ii-1)+2) = dNdx(ii,2)
                Eu(3,2*(ii-1)+1) = dNdx(ii,2)
                Eu(4,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            do ii = 1,6
                Bu(1,2*(ii-1)+1) = dNdx(ii,1)
                Bu(2,2*(ii-1)+2) = dNdx(ii,2)
                Bu(3,2*(ii-1)+1) = dNdx(ii,2)
                Bu(3,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            ! total strain variation
            delta_epsi(1:3,1) = matmul(Bu,due(uix,1))           ! Eq. 22
            
            epsi(1:3,1) = epsi_prv(1:3,1) + delta_epsi(1:3,1)   ! Eq. 22
            
            delta_e(1:4,1) = matmul(Eu,due(uix,1))              ! Eq. 22
            
            ! Eq. 1 velocity gradient
            omega(1,1:2) = delta_e(1:2,1) - delta_epsi(1:2,1)
            omega(1,3)   = delta_e(3,1) - 0.5d0*delta_epsi(3,1)
            omega(1,4)   = delta_e(4,1) - 0.5d0*delta_epsi(3,1)
            
            ! reshape
            OMEGA1(1,1) = omega(1,1)
            OMEGA1(2,2) = omega(1,2)
            OMEGA1(1,2) = omega(1,3)
            OMEGA1(2,1) = omega(1,4)
            
            ! reshape
            Ee(1,1) = delta_e(1,1)
            Ee(2,2) = delta_e(2,1)
            Ee(1,2) = delta_e(3,1)
            Ee(2,1) = delta_e(4,1)
            
            ! compute plastic strain and derivatives
            call gdndxi(ipole, xi(iintp, :),nnode/2, dNedxi)
            
            dNedx = matmul(dNedxi,invJ)
            
            Be = transpose(dNedx)
                
            call gn(ipole, xi(iintp, :),Ne)
            NeT = transpose(Ne)
            
            ! interpolation of the nodal value of the plastic strain rate variation
            delta_depsipl(1:1,1) = matmul(Ne,due(eix,1))
            
            delta_ddepsipldx(1:2,1) = matmul(Be,due(eix,1))     ! gradient of the interpolated plastic strain rate variation
            
            ! between Eqs. 17 & 18
            ! updating plastic strain from previous plastic strain rate and the time increment (IP values)
            delta_epsipl = depsipl*dt
            
            epsipl = epsipl_prv + delta_epsipl
            veps_pl = veps_pl + transpose(m)*delta_epsipl
            
            ! Hookes matrix plane strain
            H(1,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[1.0d0-nu,  nu,        0.0d0,                  nu      ]
            H(2,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        1.0d0-nu,  0.0d0,                  nu      ]
            H(3,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[0.0d0,     0.0d0,     (1.0d0-2.0d0*nu)/2.0d0, 0.0d0   ]
            H(4,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        nu,        0.0d0,                  1.0d0-nu]
            
            theta1(1,1) =  cos(twan*pi/180.d0)
            theta1(1,2) =  sin(twan*pi/180.d0)
            theta1(2,1) = -sin(twan*pi/180.d0)
            theta1(2,2) =  cos(twan*pi/180.d0)
                                    
            ag(1,2) = incriT*gamma/2.d0
            ag(2,1) = incriT*gamma/2.d0
            ag(2,2) = incriT*alphav
            
            theps = matmul(matmul(theta1,ag),transpose(theta1))
            
            epsthv(1,1) = theps(1,1)
            epsthv(2,1) = theps(2,2)
            epsthv(3,1) = theps(1,2)*2.d0
            
            dag(1,2) = DincriT*gamma/2.d0
            dag(2,1) = DincriT*gamma/2.d0
            dag(2,2) = DincriT*alphav
            
            delta_epsth = matmul(matmul(theta1,dag),transpose(theta1))
            
            delta_epsthv(1,1) = delta_epsth(1,1)
            
            delta_depsthv = delta_epsthv*dt
            
            ! adjust size and give 0 out-of-plane variation of the total strain
            delta_epsi(4,1) = 0.d0
            
            ! calculate trial stress sigma_t (Eq. 15 in Niordson 2006)
            delta_sigma_j = matmul(H,(delta_epsi-delta_epsipl*transpose(m)-delta_epsthv))
            
            ! stress and higher order stress increments
            DS = matmul(OMEGA1,SIGMA1) + transpose(matmul(OMEGA1,SIGMA1))       ! Eq.29
            
            ds1(1,1) = DS(1,1)
            ds1(1,2) = DS(2,2)
            ds1(1,3) = DS(1,2)
            ds1(1,4) = 0.d0
            
            ! classical Cauchy stress
            delta_sigma = delta_sigma_j + transpose(ds1) - sigma_prv * (delta_epsi(1,1)+delta_epsi(2,1))
            
            sigma = sigma_prv + delta_sigma
            
            ! deviatoric stress
            sigmam = (sigma(1,1)+sigma(2,1)+sigma(4,1))/3.d0
            
            dev1(1,1) = sigma(1,1) - sigmam
            dev1(2,1) = sigma(2,1) - sigmam
            dev1(3,1) = sigma(3,1)
            dev1(4,1) = sigma(4,1) - sigmam
            
            dev(1:4,1) = [dev1(1,1), dev1(2,1), dev1(3,1)*2.d0, dev1(4,1)]
            devs(1:4,1) = [dev(1,1), dev(2,1), sqrt(2.d0)*dev(3,1)/2.d0, dev(4,1)]
            
            ! von Mises stress
            sigmae = sqrt(3.d0/2.d0 *dot_product(devs(1:4,1),devs(1:4,1)))
            
            ! Elastic
            if (ssl_flag .eq. 0) then
                ! IF initialization Or first elastic step
                Ep = Ep_prv
                epsipl = epsipl_prv
                gEp = Y       ! yield strength = initial yield strength
                dg_dEp = E/n1 ! thi value is updated later; could be any value
                
                sigmac = sigmae
                Q = sigmae
                rho(1:2,1) = 0.d0
                tau(1:2,1) = 0.d0
                
                if (sigmac .eq. 0.d0) then
                    ! initialize variables before a first increment is done
                    !~call history for the first time
                    dEp = 0.d0
                    depsipl = 0.d0
                    ddepsipldx(1:2,1) = 0.d0
                    m(1,1:4) = 0.d0
                else
                    ! first elastic step
                    dEp = deps0 * (sigmac/gEp)**(1.d0/mvp)      ! Eq. 13
            
                    ! control to avoid numerical issues due to the precision of the machine/MATLAB
                    if (dEp .lt. tole) then
                        dEp = tole      ! rate of Ep dEp/dt, dt being small
                                        ! as a function of m in Eq. 13 it can be decreased
                                        ! big m -> small tole
                    end if
                    !   Updating depsipl and ddepsipl (eq. 10 and 11 in the paper)
                    ! !!!!!!! change Q to q, different notation than in the
                    ! paper (Q here corresponds to q in Borg 2006)
                    depsipl = Q * dEp / sigmac
                    ddepsipldx = (dEp/sigmac/l_int**2)*rho
                    m = 3.d0/2.d0 * transpose(dev)/sigmae
                    ssl_flag = 1        ! after first elastic step start viscoplastic case
                end if
            else
            
            ! viscoplastic computation (between Eqs.17 and 18)                
            
                delta_dEp = (depsipl/dEp) * delta_depsipl + (l_int**2/dEp) &
                            & * matmul(transpose(ddepsipldx),delta_ddepsipldx)
                    
                ! compute generalised effective stress variation (Eq. 16)
                delta_Q = (sigmac/dEp)*((mvp-1.d0)*depsipl*delta_dEp/dEp + delta_depsipl) &
                        & + ((dEp/deps0)**mvp)*dg_dEp*depsipl*dt
                        
                ! compute higher order stress variation (Eq.17)
                delta_rho_c(1:2,1) = l_int**2 *((sigmac/dEp)*((mvp-1.d0)*ddepsipldx(1:2,1)*delta_dEp(1,1)/dEp &
                            & + delta_ddepsipldx(1:2,1)) + ((dEp/deps0)**mvp)*dg_dEp*ddepsipldx(1:2,1)*dt)
                            
!               if ((delta_Q(1,1) .lt. 0.d0) .and. (Q_prv .lt. abs(delta_Q(1,1)))) then
!                      ! adjust variables
!                      ad1 = (sigmac*depsipl**2*(mvp-1)/dEp**3) + sigmac/dEp;
!                      
!                      xy =  matmul(transpose(ddepsipldx),delta_ddepsipldx)
!                     
!                      delta_depsipl_ad(1,1) = (Q_prv - (sigmac*depsipl*l_int**2 * (mvp-1.d0)/dEp**3)*xy(1,1) &
!                                         & - ((dEp/deps0)**mvp)*dg_dEp*depsipl*dt)/ad1
!                      
!                      delta_depsipl = delta_depsipl_ad;
!                      
!                      delta_dEp = (depsipl/dEp) * delta_depsipl + &
!                             & (l_int**2/dEp) * matmul(transpose(ddepsipldx),delta_ddepsipldx)
!                      
!                      ! Compute generalised effective stress variation (eq. 16 in the apaper)
!                      delta_Q = -((sigmac/dEp)*((mvp-1)*depsipl*delta_dEp/dEp + delta_depsipl) + & 
!                                & (dEp/deps0)**mvp*dg_dEp*depsipl*dt)
!                                
!                                
!                      delta_rho_c(1:2,1) = l_int**2 *((sigmac/dEp)*((mvp-1.d0)*ddepsipldx(1:2,1)*delta_dEp(1,1)/dEp &
!                           & + delta_ddepsipldx(1:2,1)) + ((dEp/deps0)**mvp)*dg_dEp*ddepsipldx(1:2,1)*dt)
!                                
!                      call mexPrintf('\n')
!                      call mexPrintf('fix active!')
!                      call mexPrintf('\n')
!               end if
                
                ! RHS of Eq. 17
                delta_rho = delta_rho_c + matmul(transpose(Ee),rho_prv)
                
                ! Eq.30
                delta_tau = delta_rho_c + matmul(transpose(Ee),tau_prv) - tau_prv*(delta_epsi(1,1)+delta_epsi(2,1))
                
                Q = Q_prv + delta_Q(1,1)
                
                rho = rho_prv + delta_rho
                
                tau = tau_prv + delta_tau
                
                ! update sigmac (Eq.12)
                sigmac = sqrt(Q**2 + dot_product(rho(1:2,1),rho(1:2,1))/l_int**2)
                
                ! update Ep_rate (Eqs. 14 & 13)
                dEp = deps0 * (sigmac/gEp_prv)**(1.d0/mvp)
                
                unloading_flag = 0
                
!!                 --- Evaluate individual conditions ---
!                 cond1 = (sigmac .lt. sigmac_prv)
!                 cond2 = (norm2(ddepsipldx) .le. deps0*1.5d0)
!                 cond3 = (dEp .lt. deps0*1.5d0)
!!                 cond2 = (norm2(ddepsipldx)*dt .le. 1e-6)
!!                 cond3 = (dEp*dt .lt. 1e-6)
!                 cond4 = (dEp*dt .lt. dEp_mem*dt)

                tol_rel_geom = 1e-7     ! relative tol
                tol_abs_geom = 1e-15    ! absolute floor in micrometres (set small)
                tol_rel_strain = 1e-7   ! relative tol 1e-6
                tol_abs_strain = 1e-15

                ! choose meaningful references
                ref1 = abs(sigmac_prv)          ! reference for sigmac comparison
                ref4 = abs(dEp_mem)             ! reference for dEp comparison

                ! compute hybrid tolerances
                tol1 = max(tol_abs_geom, tol_rel_geom * ref1)
                tol4 = max(tol_abs_strain, tol_rel_strain * ref4)

                ! Evaluate conditions with tolerance margins
                cond1 = (sigmac .lt. sigmac_prv - tol1)
                cond2 = (norm2(ddepsipldx) .le. deps0/15.d0)
                cond3 = (dEp .lt.  deps0)
                cond4 = (dEp .lt. dEp_mem - tol4)

                ! --- Set flags numerically ---
                flag_cond1 = merge(1, 0, cond1)
                flag_cond2 = merge(1, 0, cond2)
                flag_cond3 = merge(1, 0, cond3)
                flag_cond4 = merge(1, 0, cond4)

                ! --- Composite condition ---
                if (cond1 .and. cond3 .and. cond4) then
                    rho = rho_prv
                    tau = tau_prv

                    ! update sigmac (Eq.12)
                    sigmac = sqrt(Q**2 + dot_product(rho(1:2,1),rho(1:2,1))/l_int**2)

                    ! update Ep_rate (Eqs. 14 & 13)
                    dEp = deps0 * (sigmac/gEp_prv)**(1.d0/mvp)

                    unloading_flag = 1

                    call mexPrintf('\n')
                    call mexPrintf('fix active !!!!!!!!!!!!!!!!!')
                    call mexPrintf('\n')
                    
                end if
                
!                ! debug print (helpful while tuning)
!                call mexPrintf('sigmac= sigmac_prv= tol1=')
!                write(txt,'(12e15.6)') sigmac, sigmac_prv, tol1
!                call mexPrintf(txt)
!                call mexPrintf('\n')
!                
!                call mexPrintf('dEt= dEp_mem= tol4=')
!                write(txt,'(12e15.6)') dEp*dt, dEp_mem*dt, tol4
!                call mexPrintf(txt)
!                call mexPrintf('\n')
                
!!                 if ((abs(ddepsipldx(1,1)) .le. 1e-14) .and. ((dEp .lt. dEp_mem) .or. (dEp .lt. 1e-10))) then
!                if ((sigmac .lt. sigmac_prv) .and. (norm2(ddepsipldx)*dt .le. 1e-6) .and. &
!                                & (dEp*dt .lt. 1e-6) .and. (dEp*dt .lt. dEp_mem*dt)) then
!
!                    rho = rho_prv
!                    tau = tau_prv
!                    
!                    ! update sigmac (Eq.12)
!                    sigmac = sqrt(Q**2 + dot_product(rho(1:2,1),rho(1:2,1))/l_int**2)
!                    
!                    ! update Ep_rate (Eqs. 14 & 13)
!                    dEp = deps0 * (sigmac/gEp_prv)**(1.d0/mvp)
!                    
!                    unloading_flag = 1
!                    
!                    call mexPrintf('\n')
!                    call mexPrintf('fix active !!!!!!!!!!!!!!!!!')
!                    call mexPrintf('\n')
!                    
!                end if
                
                if ((sigmac/gEp_prv) .lt.0.8d0) then
                    plastic_flag = 0;
                else
                    plastic_flag = 1;
                end if
                
                ! control to avoid numerical issues due to the precision of the machine/MATLAB
                if (dEp .lt. tole) then
                    dEp = tole      ! rate of Ep dEp/dt, dt being small
                                    ! as a function of m in Eq. 13 it can be decreased
                                    ! big m -> small tole
                end if
                
                delta_Ep = dEp*dt
                
                ! update Ep
                Ep = Ep_prv + delta_Ep
                
                if (Ep .lt. 0.d0) then
                    call mexPrintf('\n')
                    call mexErrMsgTxt('Error ... Ep is negative \n')
                    call mexPrintf('\n')
                endif
                
                ! compute new hardening and stiffness
                call ghardening(Ep, matre, hardening_law, incriT, gEp, dg_dEp)
                
                ! update depsipl and ddepsipl (Eqs. 10 & 11)
                depsipl = Q * dEp / sigmac
                
                ddepsipldx = (dEp / sigmac / l_int**2) * rho
                
                ! direction of the plastic flow
                m = 3.d0/2.d0 * transpose(dev)/sigmae
                
!                 if (depsipl .lt. tole) then
!                     depsipl = tole
!                 endif
                    
                
            end if
            
            ! store history data in current
            Tempeldt(MAXNHISTI*(iintp-1)+ 13:15) = epsi(1:3,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 16)    = epsipl
            Tempeldt(MAXNHISTI*(iintp-1)+ 17:20) = sigma(1:4,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 21)    = Ep
            Tempeldt(MAXNHISTI*(iintp-1)+ 22)    = Q
            Tempeldt(MAXNHISTI*(iintp-1)+ 23:24) = tau(1:2,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 25)    = depsipl
            Tempeldt(MAXNHISTI*(iintp-1)+ 26:27) = ddepsipldx(1:2,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 28)    = dg_dEp
            Tempeldt(MAXNHISTI*(iintp-1)+ 29:32) = m(1,1:4)
            Tempeldt(MAXNHISTI*(iintp-1)+ 33)    = plastic_flag
            Tempeldt(MAXNHISTI*(iintp-1)+ 34)    = ssl_flag
            Tempeldt(MAXNHISTI*(iintp-1)+ 35)    = sigmac
            Tempeldt(MAXNHISTI*(iintp-1)+ 36)    = sigmae
            Tempeldt(MAXNHISTI*(iintp-1)+ 37)    = dEp
            Tempeldt(MAXNHISTI*(iintp-1)+ 38)    = gEp
            Tempeldt(MAXNHISTI*(iintp-1)+ 39:40) = rho(1:2,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 41:44) = veps_pl(1:4,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 45:48) = epsthv(1:4,1)
            Tempeldt(MAXNHISTI*(iintp-1)+ 49)    = dEp_mem
            Tempeldt(MAXNHISTI*(iintp-1)+ 50)    = unloading_flag
            Tempeldt(MAXNHISTI*(iintp-1)+ 51)    = flag_cond1
            Tempeldt(MAXNHISTI*(iintp-1)+ 52)    = flag_cond2
            Tempeldt(MAXNHISTI*(iintp-1)+ 53)    = flag_cond3
            Tempeldt(MAXNHISTI*(iintp-1)+ 54)    = flag_cond4
            
!             call mexPrintf('\n')
!             call mexPrintf('iintp = ')
!             write(txt,'(I2)') iintp
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('epsi = ')
!             write(txt,'(12e15.6)') epsi(1,1), epsi(2,1), epsi(3,1), epsi(4,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('epsipl =')
!             write(txt,'(12e15.6)') epsipl
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('veps_pl =')
!             write(txt,'(12e15.6)') veps_pl(1,1), veps_pl(2,1), veps_pl(3,1), veps_pl(4,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('sigma =')
!             write(txt,'(12e15.6)') sigma(1,1), sigma(2,1), sigma(3,1), sigma(4,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('Ep =')
!             write(txt,'(12e15.6)') Ep
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('Q =')
!             write(txt,'(12e15.6)') Q
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('tau =')
!             write(txt,'(12e15.6)') tau(1,1), tau(2,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('depsipl =')
!             write(txt,'(12e15.6)') depsipl
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('ddepsipldx =')
!             write(txt,'(12e15.6)') ddepsipldx(1,1), ddepsipldx(2,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('dg_dEp =')
!             write(txt,'(12e15.6)') dg_dEp
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('m =')
!             write(txt,'(12e15.6)') m(1,1), m(1,2), m(1,3), m(1,4)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('sigmac =')
!             write(txt,'(12e15.6)') sigmac
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('sigmae =')
!             write(txt,'(12e15.6)') sigmae
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('dEp =')
!             write(txt,'(12e15.6)') dEp
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('gEp =')
!             write(txt,'(12e15.6)') gEp
!             call mexPrintf(txt)
!             call mexPrintf('\n')
!             
!             call mexPrintf('rho =')
!             write(txt,'(12e15.6)') rho(1,1), rho(2,1)
!             call mexPrintf(txt)
!             call mexPrintf('\n')
            
            ! construct right hand side
            finte(1:15,1) = 0.d0
            finte(uix,1) = fu(1:12,1)
            finte(eix,1) = fe(1:3,1)
            
        end do
        
        !   output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(ndof,1))
        out1(1:ndof,1) = finte(1:ndof,1)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,(MAXNHISTI*nintp)))
        out2(1,:) = Tempeldt(1:MAXNHISTI*nintp)
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0.d0
        
        deallocate(fu)
        deallocate(fe)
        deallocate(Tempeldt)
        
        ! Computing additional term in the expression of the force (Eqs.27-28)
    elseif (key .EQ. 'fvp') then        ! viscoplastic force increment
    
        ! thickness
        t = geome
        ! Youngs Modulus
        E = matre(1,1)
        ! Elastic shear modulus
        nu = matre(1,2)
        
        if (incriT .eq. 0.d0) then
        ! rate sensitivity parameter
            mvp = matre(1,4)
        ! Intrensic length
            l_int = matre(1,3)
        else
            mvp = matre(1,15)
        ! Intrensic length
            l_int = matre(1,16)
        end if
        
        ! reference strain rate
        deps0 = matre(1,5)
        ! number of integration points
        nintp = datae(1)
        
        ! local coordinates and weight factors of intergration points
        call gintp(ipolu, nintp, xi, w)
        
        ! compute element right hand side in integration point loop
        fvpu(1:ndofu,1) = 0.d0
        fvpe(1:ndofe,1) = 0.d0
        fvp(1:ndof,1) = 0.d0

        allocate(Tempeldt(MAXNHISTI*nintp))
        Tempeldt = eldt
        
        ! initializations
        Bu(1:3,1:ndofu) = 0.d0
        H(1:4,1:4) = 0.d0
        
        do iintp = 1,nintp
        
            depsipl           = Tempeldt(MAXNHISTI*(iintp-1)+ 25)     ! rate of sclar plastic strain rate
            ddepsipldx(1:2,1) = Tempeldt(MAXNHISTI*(iintp-1)+ 26:27)  ! gradient of plastic strain rate
            dg_dEp            = Tempeldt(MAXNHISTI*(iintp-1)+ 28)     ! hardening stiffness
            m(1,1:4)          = Tempeldt(MAXNHISTI*(iintp-1)+ 29:32)  ! plastics strain direction
            dEp               = Tempeldt(MAXNHISTI*(iintp-1)+ 37)     ! rate of generalised plastic strain
            
            dEp_mem = Tempeldt(MAXNHISTI*(iintp-1)+ 49)
            
            unloading_flag = Tempeldt(MAXNHISTI*(iintp-1)+ 50)
            
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            ! strain displacement matrix Bu
            dNdx = matmul(dNdxi,invJ)
            
            do ii = 1,6
                Bu(1,2*(ii-1)+1) = dNdx(ii,1)
                Bu(2,2*(ii-1)+2) = dNdx(ii,2)
                Bu(3,2*(ii-1)+1) = dNdx(ii,2)
                Bu(3,2*(ii-1)+2) = dNdx(ii,1)
            end do
            
            ! compute plastic strain and derivatives
            call gdndxi(ipole, xi(iintp, :),nnode/2, dNedxi)
            
            dNedx = matmul(dNedxi,invJ)
            
            Be = transpose(dNedx)
                
            call gn(ipole, xi(iintp, :),Ne)
            NeT = transpose(Ne)
            
            ! Hookean matrix plane strain
            H(1,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[1.0d0-nu,  nu,        0.0d0,                  nu      ]
            H(2,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        1.0d0-nu,  0.0d0,                  nu      ]
            H(3,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[0.0d0,     0.0d0,     (1.0d0-2.0d0*nu)/2.0d0, 0.0d0   ]
            H(4,:) = E/((1.0d0+nu)*(1.0d0-2.0d0*nu))*[nu,        nu,        0.0d0,                  1.0d0-nu]
            
            ! contribution of integration points (Eq.27)
            fvpu = fvpu + dt * w(iintp) * t * depsipl * matmul(matmul(matmul(transpose(Bu),I4),H),transpose(m))*detJ
            
            qs = matmul(m,matmul(H,transpose(m)))
            
            fvpe = fvpe - dt * (w(iintp) * t * (qs(1,1) * depsipl + depsipl * dg_dEp * (dEp/deps0)**(mvp) ) * NeT * detJ &
                            & + w(iintp) * t * l_int**2 * dg_dEp * (dEp/deps0)**(mvp) * matmul(transpose(Be),ddepsipldx)*detJ)
            
            ! --- Composite condition ---
            if (unloading_flag .eq. 1) then
                
                fvpu = fvpu + dt * w(iintp) * t * deps0 * matmul(matmul(matmul(transpose(Bu),I4),H),transpose(m))*detJ
                fvpe(1:3,1) = 0.d0

                call mexPrintf('\n')
                call mexPrintf('fix fvp !!!')
                call mexPrintf('\n')
                
            end if
                            
        end do
        fvp(uix,1) = fvpu(1:12,1)
        fvp(eix,1) = fvpe(1:3,1)
        
        !   output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(ndof,1))
        out1(1:ndof,1) = fvp(1:ndof,1)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0.d0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0.d0
        
        deallocate(Tempeldt)

    !=================================================================!
    ! output control                                                  !
    !=================================================================!
    
    ! displacement
    elseif (key .eq. 'displa') then
        
        if (allocated(out1)) deallocate(out1)        
        allocate(out1(nnode,2))
        out1 = transpose(reshape(ue(uix,1), (/2,nnode/)))
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    ! plastic nodal variable Eq.23 (additional dof)
    elseif (key .eq. 'nodplv') then
        
        allocate(outn(3,1))
        outn(1:3,1) = due(eix,1)
        
        ! interpolate midside nodes
        if (allocated(out1)) deallocate(out1)
        allocate(out1(6,1))
        out1(1,1) = outn(1,1)
        out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
        out1(3,1) = outn(2,1)
        out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
        out1(5,1) = outn(3,1)
        out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
        
        deallocate(outn)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif ((key .eq. 'strain') .or. (key .eq. 'stress') .or. (key .eq. 'strapl')) then
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            if (key .eq. 'strain') then
                allocate(outi(1,4))
                outi(1,1:3) = eldt(1:3)
                outi(1,4) = 0.d0
                
            elseif (key .eq. 'stress') then
                allocate(outi(1,4))
                outi(1,1:4) = eldt(5:8)
            else
                allocate(outi(1,4))
                outi(1,1:4) = eldt(41:44)
            end if
            
            allocate (out1(6,4))
            do ii = 1,nnode
                out1(ii,1:4) = outi(1,1:4)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,4))
            fe1(1:3,1:4) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (key .eq. 'strain') then
                    allocate(outi(1,4))
                    outi(1,1:3) = eldt(MAXNHISTI*(iintp-1)+ 1:3)
                    outi(1,4) = 0.d0
                elseif (key .eq. 'stress') then
                    allocate(outi(1,4))
                    outi(1,1:4) = eldt(MAXNHISTI*(iintp-1)+ 5:8)
                else
                    allocate(outi(1,4))
                    outi(1,1:4) = eldt(MAXNHISTI*(iintp-1)+ 41:44)
                endif
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,4))
            
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,4))
            out1(1,1:4) = outn(1,1:4)
            out1(2,1:4) = (outn(1,1:4) + outn(2,1:4))/2.d0
            out1(3,1:4) = outn(2,1:4)
            out1(4,1:4) = (outn(2,1:4) + outn(3,1:4))/2.d0
            out1(5,1:4) = outn(3,1:4)
            out1(6,1:4) = (outn(3,1:4) + outn(1,1:4))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
    
    elseif (key .eq. 'rateEp') then
        
        nintp = datae(1)
        
        if (nintp .eq. 1) then
            allocate(outi(1,1))
            outi(1,1) = eldt(37)
        
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                allocate(outi(1,1))
                outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 37)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
    
    elseif ((key .eq. 'qeffec') .or. (key .eq. 'epsplr') .or. (key .eq. 'Epeffe')) then
    
        ! qeffec - generalized effective stress Eq. 6 (small q)
        ! epsplr - classical plastic strain rate Eq. 10
        ! Epeffr - generalized plastic strain Eq. 13
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            if (key .eq. 'qeffec') then
                allocate(outi(1,1))
                outi(1,1) = eldt(10)
                
            elseif (key .eq. 'epsplr') then
                allocate(outi(1,1))
                outi(1,1) = eldt(25)
            else
                allocate(outi(1,1))
                outi(1,1) = eldt(9)
            end if
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (key .eq. 'qeffec') then
                    allocate(outi(1,1))
                    outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 10)
                elseif (key .eq. 'epsplr') then
                    allocate(outi(1,1))
                    outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 25)
                else
                    allocate(outi(1,1))
                    outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 9)
                endif
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
    elseif ((key .eq. 'greplr') .or. (key .eq. 'rhohos')) then
        ! gradient of plastic strain rate and higher order stress
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            if (key .eq. 'greplr') then
                allocate(outi(1,2))
                outi(1,1:2) = eldt(26:27)
            else
                allocate(outi(1,2))
                outi(1,1:2) = eldt(39:40)
            end if
            
            allocate (out1(6,2))
            do ii = 1,nnode
                out1(ii,1:2) = outi(1,1:2)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,2))
            fe1(1:3,1:2) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (key .eq. 'greplr') then
                    allocate(outi(1,2))
                    outi(1,1:2) = eldt(MAXNHISTI*(iintp-1)+ 26:27)
                else
                    allocate(outi(1,2))
                    outi(1,1:2) = eldt(MAXNHISTI*(iintp-1)+ 39:40)
                endif
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,2))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,2))
            out1(1,1:2) = outn(1,1:2)
            out1(2,1:2) = (outn(1,1:2) + outn(2,1:2))/2.d0
            out1(3,1:2) = outn(2,1:2)
            out1(4,1:2) = (outn(2,1:2) + outn(3,1:2))/2.d0
            out1(5,1:2) = outn(3,1:2)
            out1(6,1:2) = (outn(3,1:2) + outn(1,1:2))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
    elseif ((key .eq. 'vmises') .or. (key .eq. 'sigmac')) then
        ! von Mises stress and effective stress
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            if (key .eq. 'vmises') then
                allocate(outi(1,1))
                outi(1,1) = eldt(36)
            else
                allocate(outi(1,1))
                outi(1,1) = eldt(35)
            end if
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (key .eq. 'vmises') then
                    allocate(outi(1,1))
                    outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 36)
                else
                    allocate(outi(1,1))
                    outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 35)
                endif
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
    elseif (key .eq. 'volTwn') then
        ! volume of finite element
        volume = 0.d0
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            volume = volume + w(iintp) * t * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'strnTw') then
    
        vastrn(1,1:3) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrn(1,1:3) = vastrn(1,1:3) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,3))
        out1(1,1:3) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'strsTw') then
    
        vastrs(1,1:4) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrs(1,1:4) = vastrs(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'volume') then
        ! volume of finite element
        volume = 0.d0
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            volume = volume + w(iintp) * t * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = volume
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'vastrn') then
    
        vastrn(1,1:3) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrn(1,1:3) = vastrn(1,1:3) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,3))
        out1(1,1:3) = vastrn(1,1:3)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'vaplst') then
    
        vastpl = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastpl = vastpl + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 16) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'vastrs') then
    
        vastrs(1,1:4) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrs(1,1:4) = vastrs(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = vastrs(1,1:4)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'vasigV') then
    
        sigV = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            sigV = sigV + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 36) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'vepsth') then
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        if (nintp .eq. 1) then      ! one integration point
            
            allocate(outi(1,4))
            outi(1,1:4) = eldt(45:48)
            
            allocate (out1(6,4))
            do ii = 1,nnode
                out1(ii,1:4) = outi(1,1:4)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,4))
            fe1(1:3,1:4) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                allocate(outi(1,4))
                outi(1,1:4) = eldt(MAXNHISTI*(iintp-1)+ 45:48)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,4))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,4))
            out1(1,1:4) = outn(1,1:4)
            out1(2,1:4) = (outn(1,1:4) + outn(2,1:4))/2.d0
            out1(3,1:4) = outn(2,1:4)
            out1(4,1:4) = (outn(2,1:4) + outn(3,1:4))/2.d0
            out1(5,1:4) = outn(3,1:4)
            out1(6,1:4) = (outn(3,1:4) + outn(1,1:4))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
        elseif ((key .eq. 'sstrpl')) then
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            
            allocate(outi(1,1))
            outi(1,1) = eldt(16)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                allocate(outi(1,1))
                outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 16)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
    elseif (key .eq. 'plflag') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 33)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
    elseif (key .eq. 'unflag') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 50)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
    elseif (key .eq. 'flagC1') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 51)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
    elseif (key .eq. 'flagC2') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 52)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
    elseif (key .eq. 'flagC3') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 53)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
    elseif (key .eq. 'flagC4') then
        
        nintp = datae(1)
        
        ! one integration points
        if (nintp .eq. 1) then
        
            allocate(outi(1,1))
            outi(1,1) = eldt(33)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                if (allocated (outi)) deallocate (outi)
                allocate(outi(1,1))
                outi = eldt(MAXNHISTI*(iintp-1) + 54)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
                
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        endif
        
        elseif (key .eq. 'curgep') then
    
        nintp = datae(1)
        
        if (nintp .eq. 1) then      ! one integration point
            allocate(outi(1,1))
            outi(1,1) = eldt(38)
            
            allocate (out1(6,1))
            do ii = 1,nnode
                out1(ii,1) = outi(1,1)
            end do
            
            deallocate(outi)
            deallocate(out1)
        else
            ! cross section
            t = geome
            
            ! local coordinates and weight factors of the integration points
            call gintp(ipolu, nintp, xi, w)
            
            ! build extrapolation system in integration points
            allocate (Ee1(3,3))
            Ee1(1:3,1:3) = 0.d0
            
            allocate(fe1(3,1))
            fe1(1:3,1) = 0.d0
            
            do iintp = 1,nintp
                
                ! Jacobian
                call gdndxi('t06', xi(iintp, :),nnode, dNdxi)
                J = matmul(transpose(nodee(:, 1:2)), dNdxi)
                detJ = abs(FindDet(J,2)/2.0d0)
                
                ! interpolation matrix
                call gn('t03', xi(iintp,:),Ne)
                
                !output them
                allocate(outi(1,1))
                outi(1,1) = eldt(MAXNHISTI*(iintp-1)+ 38)
                
                ! account for integration point coordinates
                Ee1 = Ee1 + w(iintp) * t * matmul(transpose(Ne),Ne) * detJ;
                fe1 = fe1 + w(iintp) * t * matmul(transpose(Ne),outi) * detJ
                deallocate(outi)
            end do
!             
            ! perform extrapolation
            allocate(outn(3,1))
            outn = matmul(inv(Ee1),fe1)
            
            deallocate(Ee1)
            deallocate(fe1)
            
            if (allocated(out1)) deallocate(out1)
            allocate(out1(6,1))
            out1(1,1) = outn(1,1)
            out1(2,1) = (outn(1,1) + outn(2,1))/2.d0
            out1(3,1) = outn(2,1)
            out1(4,1) = (outn(2,1) + outn(3,1))/2.d0
            out1(5,1) = outn(3,1)
            out1(6,1) = (outn(3,1) + outn(1,1))/2.d0
            
            deallocate(outn)
            
            if (allocated(out2)) deallocate(out2)
            allocate(out2(1,1))
            out2 = 0
            
            if (allocated(out3)) deallocate(out3)
            allocate(out3(1,1))
            out3 = 0
        end if
        
    elseif (key .eq. 'volMod') then
        ! volume of finite element
        volume = 0.d0
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            volume = volume + w(iintp) * t * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'strnMo') then
    
        vastrn(1,1:3) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrn(1,1:3) = vastrn(1,1:3) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,3))
        out1(1,1:3) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'strsMo') then
    
        vastrs(1,1:4) = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            vastrs(1,1:4) = vastrs(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
        
    elseif (key .eq. 'sigVMo') then
    
        sigV = 0.d0
        
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)
        
        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0
            
            sigV = sigV + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 36) * detJ
        end do
        
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,1))
        out1(1,1) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,1))
        out2 = 0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0
    
    elseif (key .eq. 'stsnGR') then
        ! average on IPs x volume of the element
        vastress(1,1:4) = 0.d0
        vastrain(1,1:3) = 0.d0
        volume = 0.d0
    
        ! number of integration points
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)

        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0

            ! Volume Average of Stress
            vastress(1,1:4) = vastress(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        
            ! Volume Average of Strain
            vastrain(1,1:3) = vastrain(1,1:3) + w(iintp) ** t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
            
            ! Volume  of the element
            volume = volume + w(iintp) * t * detJ

        end do

        !  output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = vastress(1,1:4)
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,3))
        out2(1,1:3) = vastrain(1,1:3)
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = volume
    
    elseif (key .eq. 'stsnTw') then
        ! average on IPs x volume of the element
        vastress(1,1:4) = 0.d0
        vastrain(1,1:3) = 0.d0
        volume = 0.d0
    
        ! number of integration points
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)

        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0

            ! Volume Average of Stress
            vastress(1,1:4) = vastress(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        
            ! Volume Average of Strain
            vastrain(1,1:3) = vastrain(1,1:3) + w(iintp) ** t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
            
            ! Volume  of the element
            volume = volume + w(iintp) * t * detJ

        end do

        !  output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,3))
        out2(1,1:3) = 0.d0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0.d0
    
    elseif (key .eq. 'stsnMo') then
        ! average on IPs x volume of the element
        vastress(1,1:4) = 0.d0
        vastrain(1,1:3) = 0.d0
        volume = 0.d0
    
        ! number of integration points
        nintp = datae(1)
        ! cross section
        t = geome
        ! local coordinates and weight factors of the integration points
        call gintp(ipolu, nintp, xi, w)

        do iintp = 1,nintp  
            ! Jacobian
            displa  = reshape(ue(uix,1), (/2, nnode/))
        
            call gdndxi(ipolu, xi(iintp, :),nnode, dNdxi)
            
            J = matmul(transpose(nodee(:,1:2))+ displa, dNdxi)
            
            invJ = inv(J)
            detJ=abs(FindDet(J,2))/2.0d0

            ! Volume Average of Stress
            vastress(1,1:4) = vastress(1,1:4) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 5:8) * detJ
        
            ! Volume Average of Strain
            vastrain(1,1:3) = vastrain(1,1:3) + w(iintp) * t * eldt(MAXNHISTI*(iintp-1)+ 1:3) * detJ
            
            ! Volume  of the element
            volume = volume + w(iintp) * t * detJ

        end do

        !  output
        if (allocated(out1)) deallocate(out1)
        allocate(out1(1,4))
        out1(1,1:4) = 0.d0
        
        if (allocated(out2)) deallocate(out2)
        allocate(out2(1,3))
        out2(1,1:3) = 0.d0
        
        if (allocated(out3)) deallocate(out3)
        allocate(out3(1,1))
        out3 = 0.d0
        
!     
    ! Illegal key
    ! -----------

    else

    !   return empty matrix
    !write(*,*) "Illegal number of integration points"
    !STOP

    if (allocated(out1)) deallocate(out1)
    allocate(out1(1,1))
    out1=0.d0

    if (allocated(out2)) deallocate(out2)
    allocate(out2(1,1))
    out2=0.d0

    if (allocated(out3)) deallocate(out3)
    allocate(out3(1,1))
    out3=0.d0
    
    endif
    
deallocate(Iden)
deallocate(I4)
deallocate(I3)

END SUBROUTINE plbr2t15Fct

END MODULE mod_plbr2t15Fct
