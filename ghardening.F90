MODULE mod_ghardening
! function to calculate hardening stress and its derivative 
! at the start of the hardening the tangent of the hardening curve
! approaches the elastic stiffness
IMPLICIT NONE

contains
SUBROUTINE ghardening(Ep, matre, hardening_law, incriT, gEp, dg_dEp)

    real*8, intent(in)          :: Ep, matre(1,14), incriT
    integer, intent(in)         :: hardening_law
    
    real*8, intent(out)         :: gEp, dg_dEp
    
    real*8                      :: E, Y, theta0, sigma_v, n
    real*8                      :: eps0
    
    E = matre(1,1)
    if (incriT .eq. 0.d0) then
        Y = matre(1,6)
    else
        Y = matre(1,6) + (matre(1,7) - matre(1,6))*incriT
    end if
    eps0 = Y/E
    
    if (hardening_law .eq. 1) then
        ! power law
        n = matre(1,8)
        
        gEp = Y * (1.d0 + Ep/eps0)**(1.d0/n)
        dg_dEp = (E/n) * (1.d0 + Ep/eps0)**(1.d0/n - 1.d0)
        
    elseif (hardening_law .eq. 2) then
        ! voce hardening law
        if (incriT .eq. 0.d0) then
            theta0 = matre(1,8)
            sigma_v = matre(1,10)
        else
            theta0 = matre(1,8) + (matre(1,9) - matre(1,8))*incriT
            sigma_v = matre(1,10) + (matre(1,11) - matre(1,10))*incriT
        end if
        
!         gEp = Y + (theta0/sigma_v)*(1.d0 - exp(-sigma_v*Ep))
!         dg_dEp = theta0*exp(-sigma_v*Ep)
        gEp = Y + sigma_v*(1.d0 - exp(-(theta0/sigma_v)*Ep))
        dg_dEp = theta0*exp(-(theta0/sigma_v)*Ep)
        
    end if
        
END SUBROUTINE ghardening

END MODULE mod_ghardening
