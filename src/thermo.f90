module mod_thermo
  !
  !use mod_common_mpi, only: istep
  !
  implicit none
  !
  real(8), parameter :: tmin = 0.d0
  !
  contains
  !
#ifdef MULT_COMP


real(8) function mass_fraction(pth,tmp,scal_liq,i)
  !
  use mod_param, only: m11,m12,m2,tmpl0,tc1,pc1,tc2,pc2        

    real(8), intent(in)        :: pth,tmp
    real(8), intent(in)        :: scal_liq
    integer, intent(in)        :: i
    real(8)                    :: pvap1,pvap2,eta1,eta2
    real(8), parameter         :: a1 = -10.03664,  &  ! Hexadecane                                                                                  
                                  b1 = +3.41426,   &  ! Hexadecane                                                                                 
                                  c1 = -6.8627,    &  ! Hexadecane 
                                  d1 = -4.8630        ! Hexadecane
   
    !real(8), parameter        :: a2 = -9.54470,   &  ! n-tetradecane                                                                                       
    !                              b2 = +3.06637,   &  ! n-tetradecane                                                                              
    !                              c2 = -6.0070,    &  ! n-tetradecane      
    !                              d2 = -4.5300        ! n-tetradecane

    real(8), parameter         :: a2 = -7.77404,  &  ! n-heptane
                                  b2 = +1.85614,  &  ! n-heptane
                                  c2 = -2.8298, &    ! n-heptane
                                  d2 = -3.5070       ! n-heptane

    !real(8), parameter        :: a2 = -7.77404,  &  ! n-heptane
    !                             b2 = +1.85614,  &  ! n-heptane
    !                             c2 = -2.8298, &    ! n-heptane
    !                             d2 = -3.5070       ! n-heptane   

    real(8), parameter         :: thr = 0.98d0
    real(8), parameter         :: eps = 1.0e-15
    real(8)                    :: mol_fraction1,mol_fraction2


    if(tmp.le.tmin) then

    eta1 = 1.d0-tmpl0/tc1
    eta2 = 1.d0-tmpl0/tc2

      else

    eta1 = 1.d0-tmp/tc1
    eta2 = 1.d0-tmp/tc2

    endif

    pvap1 = pc1*exp((a1*eta1**1.0+b1*eta1**1.5+c1*eta1**2.5+d1*eta1**5.0)/(1.d0-eta1))
    pvap2 = pc2*exp((a2*eta2**1.0+b2*eta2**1.5+c2*eta2**2.5+d2*eta2**5.0)/(1.d0-eta2))


    if(scal_liq.eq.0.d0) then
    
              mol_fraction1    = 0.d0
              mol_fraction2    = 0.d0     
      else

              mol_fraction1 = (scal_liq/m11)/(scal_liq/m11+(1-scal_liq)/m12)
              mol_fraction2 = 1 - mol_fraction1
              
    endif
   
   
    select case(i)  
    
    case(1)
   
    mass_fraction = (pvap1*m11*mol_fraction1)/(pth*m2+mol_fraction1*pvap1*(m11-m2)+mol_fraction2*pvap2*(m12-m2))
       
    case(2)

    mass_fraction = (pvap2*m12*mol_fraction2)/(pth*m2+mol_fraction2*pvap2*(m12-m2)+mol_fraction1*pvap1*(m11-m2))

    end select         
 
    !mass_fraction = max(eps,min(mass_fraction,thr))       
    return 

  end function mass_fraction


#else
  real(8) function mass_fraction(pth,tmp)
    !
    use mod_param, only: m11,m2,tmpl0,tc1,pc1
    !
    real(8), intent(in) :: pth,tmp
    !
    real(8), parameter :: a = -7.77404,  &  ! n-heptane
                          b = +1.85614,  &  ! n-heptane
                          c = -2.8298, &    ! n-heptane
                          d = -3.5070       ! n-heptane
    !
    real(8), parameter :: thr = 0.98d0
    real(8), parameter :: eps = 1.0e-09
    !
    real(8) :: pvap,eta
    !
    if(tmp.le.tmin) then
      eta = 1.d0-tmpl0/tc1
    else
      eta = 1.d0-tmp  /tc1
    endif
    pvap = pc1*exp((a*eta**1.0+b*eta**1.5+c*eta**2.5+d*eta**5.0)/(1.d0-eta))
    mass_fraction = pvap*m11/(pvap*m11+(pth-pvap)*m2)
    mass_fraction = max(eps,min(mass_fraction,thr))
    !
    return
  end function mass_fraction
#endif
 

  
  !
  real(8) function thermo_rhog(pth,tmp,yl)
    !
#ifdef LOW_MACH
    use mod_param, only: ru,m11,m2,tmpg0
#else
    use mod_param, only: rho2_0
#endif
    use mod_param, only: tmpg0
    !
    real(8), intent(in) :: pth,tmp,yl
    !
#ifdef LOW_MACH
    if(tmp.le.tmin) then
      thermo_rhog = ((yl/m11+(1.d0-yl)/m2)**(-1.d0))*pth/(ru*tmpg0) ! 3rd level of approx
    else 
      thermo_rhog = ((yl/m11+(1.d0-yl)/m2)**(-1.d0))*pth/(ru*tmp  ) ! 3rd level of approx
    endif
#else
    thermo_rhog = rho2_0
#endif
    !
    return
  end function thermo_rhog
  !
  real(8) function thermo_d_lg(pth,tmp,yl,i)
    !
    ! modified Wilke Lee's law
    !
    use mod_param, only: d_m12,d_m12_2,tmpg0,pth_0,ru,m11,m2,mu2
    !
    real(8), intent(in) :: pth,tmp,yl
    integer, intent(in)        :: i
    !
#ifdef VAR_D_LG
    if(tmp.le.tmin) then
      thermo_d_lg = d_m12*((tmpg0/tmpg0)**(3.d0/2.d0))*(pth_0/pth)
    else
      thermo_d_lg = d_m12*((tmp  /tmpg0)**(3.d0/2.d0))*(pth_0/pth)
    endif
#else
    select case(i)
    case(1)
    thermo_d_lg = d_m12
    case(2)
    thermo_d_lg = d_m12_2
    end select
#endif
    ! 
    return
  end function thermo_d_lg
  !
  real(8) function thermo_mug(tmp)
    !
    ! modified Sutherland's law
    ! 
    use mod_param, only: tmpg0,mu2
    !
    real(8), intent(in) :: tmp
    !
#ifdef VAR_MUG
    if(tmp.le.tmin) then
      thermo_mug = mu2*((tmpg0/tmpg0)**(2.d0/3.d0))
    else
      thermo_mug = mu2*((tmp  /tmpg0)**(2.d0/3.d0))
    endif
#else
    thermo_mug = mu2
#endif
    !
    return
  end function thermo_mug
  !
  real(8) function thermo_kag(tmp)
    !
    ! kappa from modified Sutherland's law to ensure the same Pr
    ! 
    use mod_param, only: tmpg0,mu2,cp2,kappa2
    !
    real(8), intent(in) :: tmp
    !
#ifdef VAR_KAG
    if(tmp.le.tmin) then
      thermo_kag = (mu2*((tmpg0/tmpg0)**(2.d0/3.d0)))*cp2/pr
    else
      thermo_kag = (mu2*((tmp  /tmpg0)**(2.d0/3.d0)))*cp2/pr
    endif
#else
    thermo_kag = kappa2
#endif
    !
    return
  end function thermo_kag
  ! 
end module mod_thermo
