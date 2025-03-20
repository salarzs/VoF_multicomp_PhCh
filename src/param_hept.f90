module mod_param
  !
  implicit none
  !
  include 'bc.h90'
  !
  ! parameters for the grid
  !
  integer, parameter :: ndims = 2
  real(8), parameter :: pi = acos(-1.)
  integer, dimension(ndims), parameter :: dims = (/1,1/)
#ifdef TWOD
  integer, parameter :: itot = 2, jtot = 128, ktot = 128
#else
  integer, parameter :: itot = 128, jtot = 128, ktot = 128
#endif
  integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: n  = (/imax,jmax,kmax/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
  real(8), parameter :: ly = 0.001d0
  real(8), parameter :: lz = 0.001d0
#ifdef TWOD
  real(8), parameter :: lx = (1.0*itot/jtot)*ly
#else
  real(8), parameter :: lx = 0.001d0
#endif
  !
  real(8), parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  real(8), parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
  real(8), parameter, dimension(3) :: dl  = (/dx ,dy ,dz /)
  !
  ! simulations checks
  !
  integer, parameter :: nstep = 100
  !
  integer, parameter :: istep_do_av = 1e+7
  integer, parameter :: isave       = 1000000
  !
  logical, parameter :: start_f_ze      = .true.  ! start from 0
  !logical, parameter :: start_f_ze      = .false.  ! start from 0
  !
  logical, parameter :: restart_sp      = .false.  ! re-start from a pre-initialized single phase field
  !logical, parameter :: restart_sp      = .true.  ! re-start from a pre-initialized single phase field
  !
  logical, parameter :: restart_mp_isot = .false.  ! re-start from a pre-initialized isothermal multiphase field (no ht, no mt)
  !logical, parameter :: restart_mp_isot = .true.  ! re-start from a pre-initialized isothermal multiphase field (no ht, no mt)
  !
  logical, parameter :: restart_mp_comp = .false.  ! re-start from the full multiphase (with heat and mass) transfer
  !logical, parameter :: restart_mp_comp = .true.   ! re-start from the full multiphase (with heat and mass) transfer
  !
  logical, parameter :: var_dt      = .true.
  real(8), parameter :: dt_input    = 1e-6
  !
  ! dimensionless physical parameters
  !
  !real(8), parameter :: re  = 100.d0          
  !real(8), parameter :: re  = 64.d0          
  !real(8), parameter :: re  = 24.d0          ! droplet diameter
  real(8), parameter :: re  = 1.15d0         ! droplet diameter
  real(8), parameter :: pr  = 0.71d0           
  real(8), parameter :: we  = 0.10d0           
  real(8), parameter :: sc  = 0.87d0         ! eva
  real(8), parameter :: ste = 5.00d0         ! eva! scaled by gas temperature
  real(8), parameter :: lambda_rho0 = 25.0d0  
  real(8), parameter :: lambda_mu0  = 8.00d0  
  real(8), parameter :: lambda_cp0  = 2.22d0  
  real(8), parameter :: lambda_ka0  = 3.63d0  
  real(8), parameter :: lambda_mm1  = 3.44d0 ! [molar mass of n-Heptane/molar mass of dry air],0.1 
  real(8), parameter :: lambda_mm2  = 3.44d0 ! [molar mass of n-Heptane/molar mass of dry air],0.1 
  real(8), parameter :: lambda_mav  = 3.44d0 ! [molar mass of n-Heptane/molar mass of dry air],0.1 
  !real(8), parameter :: lambda_mm1  = 7.80d0  !226
  !real(8), parameter :: lambda_mm2  = 6.83d0  !198
  !real(8), parameter :: lambda_mav  = 7.32d0  !198
  real(8), parameter :: lambda_dcp0 = 1.00d0 ! eva
  real(8), parameter :: phi_p1      = 0.20d0 ! eva (ru/(mm2*cp2))
  real(8), parameter :: phi_p2      = 1.00d0 ! eva (pth0/((ru/mm2)*tmp0))
  real(8), parameter :: phi_p3      = 1.29d0 ! eva (pth_0/pc) 
  real(8), parameter :: phi_p4      = 1.50d0 ! eva (tmp0/tc) 
  !real(8), parameter :: phi_p4      = 1.00d0 ! eva (tmp0/tc) 
  !real(8), parameter :: phi_p4      = 0.75d0 ! eva (tmp0/tc) 
  !
  ! referencie
  !
  real(8), parameter :: k0_wave = 2.d0*pi/lx
  integer, parameter :: k0_freq = 1
  real(8), parameter, dimension(3) :: gacc = (/0.d0,0.d0,0.d0/)
  real(8), parameter :: tc1 = 540.61, pc1 = 2.736e+6
  real(8), parameter :: tc2 = 540.61, pc2 = 2.736e+6
  real(8), parameter :: ru = 8314.d0
  real(8), parameter :: m2 = 28.9647d0
  !real(8), parameter :: res    = 50.d0
  real(8), parameter :: res    = 12.d0
  real(8), parameter :: lref   = 0.25d0*lz !res*dx
  real(8), parameter :: nu2    = 8.76e-07
  !real(8), parameter :: nu2    = 1.0/re
  real(8), parameter :: d_m12    = 1e-7!nu2/sc
  real(8), parameter :: d_m12_2  = 1e-7!nu2/sc
  !real(8), parameter :: uref   = d_m12/lref
  real(8), parameter :: uref   = d_m12/lref!re*k0_wave*nu2
  !real(8), parameter :: uref   = 1.0
  real(8), parameter :: u_abc  = uref
  real(8), parameter :: pth_0  = 30*101325!phi_p3*pc
  real(8), parameter :: tmp0   = 400!phi_p4*tc
  real(8), parameter :: rho2_0 = pth_0/(phi_p2*(ru/m2)*tmp0)
  real(8), parameter :: beta_l_th = 1.d0/tmp0
  real(8), parameter :: beta_g_th = 1.d0/tmp0
  !
  real(8), parameter :: tmpg0 = tmp0                                   !same tmp for both phases!!!!
  real(8), parameter :: tmpl0 = tmpg0!431.75d0  ! phi_p4 = 1.50
  !
  ! initial flow condition
  !
  character(len=3), parameter :: inivel = 'zer' 
  character(len=3), parameter :: inivof = 'mub' 
  character(len=3), parameter :: initmp = 'sin'
  character(len=3), parameter :: inisca = 'std'
  !
  ! output frequency (so we know t_ff)
  !
  integer, parameter :: icheck    = 1
  !integer, parameter :: iout0d    = 50
  integer, parameter :: iout0d    = 1
  integer, parameter :: iout0d_ta = 5000
  integer, parameter :: iout_av   = 1000    ! averaging frequency (not necesseraly the printing one)
  integer, parameter :: iout1d    = 1000*20 ! print 
  integer, parameter :: iout2d    = 1000*20
  !integer, parameter :: iout3d    = 5000*40
  integer, parameter :: iout3d    =  1
  !integer, parameter :: iout3d    = 500000
  !integer, parameter :: iout3d    = 1
  !integer, parameter :: iout3d    = 20
  !integer, parameter :: iout3d    = 10
  !
  !integer, parameter :: icheck  = 10
  !integer, parameter :: iout0d  = 50*4!*2
  !integer, parameter :: iout1d  = int(frac*t_ff/dt) ! print
  !integer, parameter :: iout_av = 100   ! averaging frequency (not necesseraly the printing one)
  !integer, parameter :: iout2d  = 1000*20!*5
  !integer, parameter :: iout3d  = 5000*40!*5
  !
  ! dimensional quantities
  !
  real(8), parameter :: m11       = lambda_mm1*m2
  real(8), parameter :: mav       = lambda_mav*m2
  real(8), parameter :: m12       = lambda_mm2*m2
  real(8), parameter :: cp2       = 1.0086e+3!ru/(phi_p1*m2)
  real(8), parameter :: mu2       = 2.303e-5!nu2*rho2_0
  real(8), parameter :: rho1      = lambda_rho0*rho2_0
  real(8), parameter :: sigmaca   = rho2_0*(uref**2.d0)*lref/we ! we use the density of the liquid like Lohse-Verzicco
  real(8), parameter :: mu1       = lambda_mu0*mu2
  real(8), parameter :: kappa2    = 0.033d0!mu2*cp2/pr
  real(8), parameter :: kappa1    = lambda_ka0*kappa2
  real(8), parameter :: cp1       = lambda_cp0*cp2
  !real(8), parameter :: d_m12     = (mu2/rho2_0)*(1.d0/sc)
  real(8), parameter :: lheat     = 3.5051e+05!cp2*tmp0/ste ! be sure to put a good Ste
  real(8), parameter :: delta_cp  = cp1-cp2!lambda_dcp0*cp2
  !
  integer, parameter :: e         = 0
  real(8), parameter :: sinit     = 0.d0
  real(8), parameter :: gri       = 0.d0
  real(8), parameter :: cfl       = 0.50d0 !different compare to previous version
  real(8), parameter :: cfl_o     = 0.90d0
  real(8), parameter :: small     = 1e-09
  real(8), parameter :: theta_thr = 0.25d0
  real(8), parameter :: mflux0    = 0.1d0!0.0d0
  !
  ! disperse phase 
  !
  integer, parameter :: n_dp  = 1 ! initial number of droplets
  real(8), parameter :: xcc(1:n_dp) = (/0.500d0/)*lx
  real(8), parameter :: ycc(1:n_dp) = (/0.500d0/)*ly
  real(8), parameter :: zcc(1:n_dp) = (/0.500d0/)*lz
  real(8), parameter :: rd(1:n_dp)  = lref/2.d0
  !logical, parameter :: manual_dp   = .false.
  logical, parameter :: manual_dp   = .true.
  !
  character(len=*), parameter :: datadir = 'data/'
  character(len=*), parameter :: datapos = 'data/post/'
  character(len=*), parameter :: datapos_mass = 'data/post/mass/'
  character(len=*), parameter :: datapos_vol  = 'data/post/vol/'
  !
  logical, parameter, dimension(2,3) :: no_outflow = & 
      reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                .false.,.false.,   & ! no outflow in y lower,upper bound
                .false.,.false./), & ! no outflow in z lower,upper bound
                 shape(no_outflow))
  !
  ! re-calculation of the dimensionless physical parameters
  !
  real(8), parameter :: re_c          = rho2_0*uref/(mu2*k0_wave)
  real(8), parameter :: pr_c          = mu2*cp2/kappa2
  real(8), parameter :: we_c          = rho2_0*(uref**2.d0)*lref/sigmaca
  real(8), parameter :: sc_c          = nu2/d_m12                     ! only eva
  real(8), parameter :: ste_c         = cp2*tmp0/lheat                ! only eva
  real(8), parameter :: lambda_rho0_c = rho1/rho2_0                  
  real(8), parameter :: lambda_mu0_c  = mu1/mu2                      
  real(8), parameter :: lambda_cp0_c  = cp1/cp2                      
  real(8), parameter :: lambda_ka0_c  = kappa1/kappa2                
  real(8), parameter :: lambda_mm0_c  = m11/m2                         ! only eva
  real(8), parameter :: lambda_dcp0_c = delta_cp/cp2                  ! only eva
  real(8), parameter :: phi_p1_c      = (ru/(m2*cp2))                 ! only eva
  real(8), parameter :: phi_p2_c      = (pth_0/(rho2_0*(ru/m2)*tmp0))       
  real(8), parameter :: phi_p3_c      = pth_0/pc1                      ! only eva
  real(8), parameter :: phi_p4_c      = tmp0/tc1                       ! only eva
  !
  ! forcing
  !
  character(len=3), parameter :: turb_type = 'abc' 
  !character(len=3), parameter :: c_or_t    = 'tdp' 
  character(len=3), parameter :: c_or_t    = 'uni' 
  real(8)         , parameter, dimension(3) :: abc = (/u_abc,u_abc,u_abc/)
  real(8)         , parameter :: f0_t = nu2*(k0_wave)**2
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 1.d0    ! m1 - n1
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 10.d0   ! m1 - n2
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 50.d0   ! m1 - n3
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 100.d0  ! m1 - n3
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 1.d0    ! m2 - n1
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 10.d0   ! m2 - n2
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 50.d0   ! m2 - n3
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 100.d0  ! m2 - n3
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 20.0d0, amp_n = 1.d0    ! m3 - n1
  real(8)         , parameter :: amp_a = 0.d0, amp_b = 20.0d0, amp_n = 1.d0    ! m3 - n1
  !
#ifdef MULT_COMP
  integer, parameter :: or_ex = 2, it_ex = int(20.d0/theta_thr)
  real(8), parameter :: diff_coef_1_12m = 4.5e-05  
#endif
  !
end module mod_param
