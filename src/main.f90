!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!         CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS  (*) !
!      CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S     !
!    CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S     !
!   C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS     !
!  C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S                 !
! C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S                 !
! C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS              !
! C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS         !
! C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS       !
! C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S      !
! C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S     !
!  C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S     !
!   C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S     !
!    CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S     !
!      CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS      !
!         CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS        !
! * Two-Fluid with low-Mach in the gas phase                                                !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
! Two-Fluid code based on                                                                   !
! CaNS -- Canonical Navier-Stokes Solver (github.com/p-costa/CaNS)                          !
! Pedro Costa (p.simoes.costa@gmail.com)                                                    !
!                                                                                           !
! Original VoF module (vof.f90) implemented by Marco Edoardo Rosti.                         !
!                                                                                           !
! The pressure-splitting method that allows for a FFT-based Poisson solver follows the      !
! following reference:                                                                      !
!                                                                                           !
!  * Frantzis, C., & Grigoriadis, D. G. E. (2019). An efficient method for two-fluid        !
!    incompressible flows appropriate for the immersed boundary method.                     !
!    Journal of Computational Physics, 376, 28-53                                           !
!                                                                                           !
! which improves the temporal accuracy of the pressure-splitting technique described in:    !
!                                                                                           !
!  * Dodd, M. S., & Ferrante, A. (2014). A fast pressure-correction method for              !
!    incompressible two-fluid flows.                                                        !
!    Journal of Computational Physics, 273, 416-434.                                        !
!                                                                                           !
! Other references are added in the modules when needed.                                    !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!   
program cans
  !
  use iso_c_binding, only: C_PTR
  use mod_param
  use decomp_2d
  use mod_common_mpi!, only: istep
  use mod_cmpt_divth
  use mod_initmpi
  use mod_initvof
  use mod_initflow
  use mod_initgrid
  use mod_bound
  use mod_chkdiv
  use mod_chkdt
  use mod_load
  use mod_rk
  use mod_rks
  use mod_robin_bc
  use mod_initsolver
  use mod_fillps      !, only: fillps,fillpslm1,fillpslm2,restpslm
#ifdef MULTI_GRID
  use mod_solver_vc
#else
  use mod_solver
  use mod_fft         , only: fftini,fftend
#endif
  use mod_correc
  use mod_output
  use mod_vof
  use mod_sanity
  use mod_thermo
  use mod_inte_vel
  use mod_vof_to_ls
  !
  implicit none
  !
  real(8) ,dimension(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3)   :: u,v,w
  real(8) ,dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: p,pold,pp
  real(8) ,dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: up,vp,wp
  real(8) ,dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: vof,d_thinc,curv
  real(8) ,dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: mu,rho,kappa!,cpp
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1,3) :: nor
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1,6) :: cur
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: dudtrkold,dvdtrkold,dwdtrkold,dtmpdtrkold
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: div_th,divg_th!,div_s,e_div_v
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: rho_gas,rho_liq,d_lg,mu_gas,kappa_gas,g_sn
  real(8), dimension(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3)   :: phi
  real(8), dimension(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3)   :: tmp
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: ug,vg,wg,mflux,mflux2,aux
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: sca,scae,tmpge,tmple, &
                                                         normlx,normly,normlz,delta!,sc_f,pr_f
#ifdef MULT_COMP
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: sca_liq,sca2_liq,sca2e_liq,scae_liq,sca2,sca2e,d_l12,d_l2g, &
                                                         sca_test, scae_test
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: j_n,j2_n,mflux_a,mflux2_a,mflux_tot,g_sn2
#endif
  !
  real(8), dimension(-2:n(3)+3) :: dzfi,dzci,dzc,dzf,asr,zc,zf
  real(8) :: pth,ptho,dpthdt_n,dpthdt
  real(8) :: mfxt,mgas, v_left,w_left,v_bottom,w_bottom,p_left,p_bottom,u_bottom 
  !
  integer, parameter        :: num    = +3
  integer, parameter        :: rdir_i = +1
  real(8), dimension(num+1) :: tinf_v,yinf_v
  !
#ifdef USE_VOF
  integer :: num_dp,ww
  real(8), allocatable, dimension(:,:) :: pos
  real(8), allocatable, dimension(:)   :: xc_p,yc_p,zc_p,rc_p
#endif
  !
  real(8), dimension( 0:n(1)+1, 0:n(2)+1, 0:n(3)+1)   :: v_1,v_2
  real(8), dimension(10) :: var
  real(8)                :: nondim
  !character(len=7) :: fldnum
  character(len=12) :: fldnum
  real(8) :: dt,dto,time,vol0
  real(8) :: e_div,e_div_mean,divu
  real(8) :: rho0_min,d_lg_max,mu_g_max,ka_g_max
  real(8) :: dtmax,dtmax_o
  integer :: istep,istep0,i_av
  !integer :: istep0,i_av
  integer :: k,l,i,j,kk,lenr
  real(8) :: divmax,divtot
  character(len=1) :: action_load
#ifdef TWOD
  real(8) :: u_max,u_min
#endif
  real(8) :: v_max,v_min,last_saved
  logical :: is_first_tmp,is_first_vel,is_fir_field,kill,is_data,is_first_pos,do_avg
  !
  integer, parameter :: qb = abs(lbound(u ,1))
  integer, parameter :: qs = abs(lbound(ug,1))
  !
#ifdef TIMING
  real(8) :: dt12,dt12av,dt12min,dt12max
#endif
  !character(len=7) :: istepchar
  !
#ifndef MULTI_GRID
  type(C_PTR), dimension(2,2) :: arrplanp
  real(8)    , dimension(imax,jmax) :: lambdaxyp
  real(8)    , dimension(ktot) :: ap,bp,cp
  real(8) :: normfftp
  type rhs_bound
    real(8), dimension(n(2),n(3),0:1) :: x
    real(8), dimension(n(1),n(3),0:1) :: y
    real(8), dimension(n(1),n(2),0:1) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
#endif
  !
  ! create directories
  !
  inquire(file='data/',exist=is_data)
  if(.not.is_data) call execute_command_line('mkdir -p data')
  inquire(file='data/post/',exist=is_data)
  if(.not.is_data) call execute_command_line('mkdir -p data/post')
  inquire(file='data/post/mass/',exist=is_data)
  if(.not.is_data) call execute_command_line('mkdir -p data/post/mass')
  inquire(file='data/post/vol/',exist=is_data)
  if(.not.is_data) call execute_command_line('mkdir -p data/post/vol')
  !
  rho_liq(:,:,:) = rho1
  !
  call initmpi(ng,cbcpre)
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, ''
  !
  ! grid generation along z 
  !
  call initgrid(inivel,n(3),gri,lz,dzc,dzf,zc,zf,asr)
  if(myid.eq.0) then
    inquire (iolength=lenr) dzc(1)
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*n(3)*lenr)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=-2,ktot+3
      write(99,'(6E15.7)') 1.d0*kk,zf(kk),zc(kk),dzf(kk),dzc(kk),asr(kk)
    enddo
    close(99)
  endif
  !
  dzci = dzc**(-1.d0)
  dzfi = dzf**(-1.d0)
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,n,dims,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                   dli,dzci,dzfi)
  !
  ! initialize Poisson/Helmholtz solver
  !
#ifndef MULTI_GRID
  call initsolver(n,dli,dzci,dzfi,cbcpre,bcpre(:,:),lambdaxyp,(/'c','c','c'/),ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
#else
  call init_solver_mg(n,dli,dzci,cbcpre,bcpre)
#endif
  !
#ifdef USE_VOF
  if(.not.manual_dp) then
    !
    ! load number and positions of droplets
    !
    if(myid.eq.0) call system('rm -f preprocessing/n_dp')
    if(myid.eq.0) call system('cat preprocessing/num_dp.out|wc -l>>preprocessing/n_dp')
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !
    open(29,file='preprocessing/n_dp',status='old',action='read',iostat=ierr)
     read(29,*) num_dp
    close(29)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !
    allocate(pos(num_dp,4),xc_p(num_dp),yc_p(num_dp),zc_p(num_dp),rc_p(num_dp))
    open(30,file='preprocessing/num_dp.out',status='old',action='read',iostat=ierr)
     do ww = 1,num_dp
       read(30,*) pos(ww,1),pos(ww,2),pos(ww,3),pos(ww,4)
     enddo
    close(30)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !
    do ww=1,num_dp
      xc_p(ww) = pos(ww,1)
      yc_p(ww) = pos(ww,2) 
      zc_p(ww) = pos(ww,3)
      rc_p(ww) = pos(ww,4) ! already radius
    enddo
    !
  else
    !
    ! do it manually
    !
    allocate(pos(n_dp,4),xc_p(n_dp),yc_p(n_dp),zc_p(n_dp),rc_p(n_dp))
    xc_p(1:n_dp) = xcc(1:n_dp)
    yc_p(1:n_dp) = ycc(1:n_dp)
    zc_p(1:n_dp) = zcc(1:n_dp)
    rc_p(1:n_dp) = rd( 1:n_dp)
    num_dp       = n_dp
    !
  endif
  !
  if(myid.eq.0) then
    print*, "initial volume:  ", num_dp*(4.d0*pi/3.d0)*minval(rc_p(:))**3
    print*, "initial # of dp: ", num_dp
  endif
#endif
  !
  if(start_f_ze) then
    !
    if(myid.eq.0) print*, "============================================="
    if(myid.eq.0) print*, " I start from 0 the simulation "
    if(myid.eq.0) print*, "============================================="
    !
    istep = 0
    time  = 0.d0
    is_first_vel = .true.
    is_first_tmp = .true.
    is_first_pos = .true.
    !
#ifdef USE_VOF
    !call initvof(n,dli,vof)
    call initvof_ma(n,dli,num_dp,xc_p,yc_p,zc_p,rc_p,vof)
#else
    vof(:,:,:) = 0.d0
#endif
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
    !
    call initvel(inivel,n,zc/lz,dzc/lz,dzf/lz,0.d0,0.00d0,time,vof,u,v,w,p)
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    !
    pold(   :,:,:) = p(:,:,:)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pold)
    !
    dudtrkold(:,:,:) = 0.d0
    dvdtrkold(:,:,:) = 0.d0
    dwdtrkold(:,:,:) = 0.d0
    !
#ifdef USE_VOF
    call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
#endif
#ifdef ENERGY
    call inittmp(initmp,n,dli,phi,vof,tmp)
#else
    tmp(:,:,:) = tmp0
#endif
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    ! 
    pth      = pth_0
    dpthdt_n = 0.d0
    dpthdt   = pth*dpthdt_n
    !
#ifdef USE_VOF
    call update_vof(n,dli,-qb,dzc,dzf,vof,nor,cur,curv,d_thinc)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          rho_gas(i,j,k)   = thermo_rhog(pth,tmp(i,j,k),sca(i,j,k))
          mu_gas(i,j,k)    = thermo_mug(tmp(i,j,k))
          kappa_gas(i,j,k) = thermo_kag(tmp(i,j,k))
          d_lg(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca(i,j,k),1)
#ifdef MULT_COMP       
          d_l2g(i,j,k)     = thermo_d_lg(pth,tmp(i,j,k),sca2(i,j,k),2)
          d_l12(i,j,k)     = diff_coef_1_12m
#endif
          rho(i,j,k)   = rho1*vof(  i,j,k)+(1.d0-vof(i,j,k))*rho_gas(  i,j,k)
          mu(i,j,k)    = mu1*vof(   i,j,k)+(1.d0-vof(i,j,k))*mu_gas(   i,j,k)
          kappa(i,j,k) = kappa1*vof(i,j,k)+(1.d0-vof(i,j,k))*kappa_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho_gas  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu_gas   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa_gas)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_lg     )
#ifdef MULT_COMP
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l2g    )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l12    )
#endif 
    !
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
#else
    rho(0:n(1)+1,0:n(2)+1,0:n(3)+1)       = rho2_0
    mu(0:n(1)+1,0:n(2)+1,0:n(3)+1)        = mu2
    kappa(0:n(1)+1,0:n(2)+1,0:n(3)+1)     = kappa2
    rho_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)   = rho2_0
    mu_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)    = mu2
    kappa_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1) = kappa2
    d_lg(0:n(1)+1,0:n(2)+1,0:n(3)+1)      = d_m12
#endif
    !
#ifdef LOW_MACH
    rho0_min = minval(rho(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,rho0_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
    rho0_min = rho2_0 ! gas phase 
#endif
    !
#ifdef VAP_MASS
    !
    !call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
    call cmpt_normals(n,dli,phi,normlx,normly,normlz)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlx)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normly)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlz)
    !
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmple,-1.d0,0.d0) ! Liq. tmp. extension
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmpge,+1.d0,0.d0) ! Gas. tmp. extension
    !
#ifdef MULT_COMP

    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          !
          if(phi(i,j,k).gt.0.d0) then           
            sca_liq(i,j,k)  = 0.90d0
          else
            sca_liq(i,j,k)  = 0.d0
          endif
        enddo
      enddo
    enddo
   
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca_liq)
    !
    call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,sca_liq,scae_liq,-1.d0,0.d0)
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca,scae_liq,1)
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_l2g, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca2,scae_liq,2)   
#else
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca)
#endif     
    !
#endif
    !
    ug(:,:,:) = 0.d0
    vg(:,:,:) = 0.d0
    wg(:,:,:) = 0.d0
    ! 
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
    !
  elseif(restart_sp) then
    !
    if(myid.eq.0) print*, "============================================="
    if(myid.eq.0) print*, " I restart by loading the single phase field "
    if(myid.eq.0) print*, "============================================="
    !
    is_first_vel = .true.
    is_first_tmp = .true.
    is_first_pos = .true.
    !
    ! load files
    !
    action_load = 'r' 
    call load(action_load,trim(datadir)//'fldu.bin'  ,n,u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldv.bin'  ,n,v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldw.bin'  ,n,w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldp.bin'  ,n,p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvof.bin',n,vof(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldug.bin' ,n,ug(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvg.bin' ,n,vg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldwg.bin' ,n,wg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldtmp.bin',n,tmp(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldsca.bin',n,sca(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(datadir)//'scalar.out',pth,dpthdt_n,time,istep)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at = ', time, 'time step = ', istep, '. ***'
    !
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    !
#ifdef USE_VOF
    !call initvof(n,dli,vof)
    call initvof_ma(n,dli,num_dp,xc_p,yc_p,zc_p,rc_p,vof)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
#endif
    !
    pold(   :,:,:) = p(:,:,:)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pold)
    !
    dudtrkold(:,:,:) = 0.d0
    dvdtrkold(:,:,:) = 0.d0
    dwdtrkold(:,:,:) = 0.d0
    !
#ifdef USE_VOF
    call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
#endif
#ifdef ENERGY
    call inittmp(initmp,n,dli,phi,vof,tmp)
#else
    tmp(:,:,:) = tmp0
#endif
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    ! 
    pth      = pth_0
    dpthdt_n = 0.d0
    dpthdt   = pth*dpthdt_n
    !
#ifdef USE_VOF
    call update_vof(n,dli,-qb,dzc,dzf,vof,nor,cur,curv,d_thinc)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          rho_gas(i,j,k)   = thermo_rhog(pth,tmp(i,j,k),sca(i,j,k))
          mu_gas(i,j,k)    = thermo_mug(tmp(i,j,k))
          kappa_gas(i,j,k) = thermo_kag(tmp(i,j,k))
          d_lg(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca(i,j,k),1)
#ifdef MULT_COMP
          d_l2g(i,j,k)     = thermo_d_lg(pth,tmp(i,j,k),sca2(i,j,k),2)
          d_l12(i,j,k)     = diff_coef_1_12m
#endif
          rho(i,j,k)   = rho1*vof(  i,j,k)+(1.d0-vof(i,j,k))*rho_gas(  i,j,k)
          mu(i,j,k)    = mu1*vof(   i,j,k)+(1.d0-vof(i,j,k))*mu_gas(   i,j,k)
          kappa(i,j,k) = kappa1*vof(i,j,k)+(1.d0-vof(i,j,k))*kappa_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho_gas  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu_gas   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa_gas)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_lg     )
    !
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
#else
    rho(0:n(1)+1,0:n(2)+1,0:n(3)+1)       = rho2_0
    mu(0:n(1)+1,0:n(2)+1,0:n(3)+1)        = mu2
    kappa(0:n(1)+1,0:n(2)+1,0:n(3)+1)     = kappa2
    rho_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)   = rho2_0
    mu_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)    = mu2
    kappa_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1) = kappa2
    d_lg(0:n(1)+1,0:n(2)+1,0:n(3)+1)      = d_m12
#endif
    !
#ifdef LOW_MACH
    rho0_min = minval(rho(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,rho0_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
    rho0_min = rho2_0 ! gas phase 
#endif
    !
#ifdef VAP_MASS
    !
    !call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
    call cmpt_normals(n,dli,phi,normlx,normly,normlz)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlx)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normly)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlz)
    !
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmple,-1.d0,0.d0) ! Liq. tmp. extension
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmpge,+1.d0,0.d0) ! Gas. tmp. extension
    !
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca &
#ifdef MULT_COMP
                  ,scae_liq,1)
#else
)
#endif
#ifdef MULT_COMP    
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca2,scae_liq,2)    !new
#endif
#endif
    !
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca)
    !
    ug(:,:,:) = 0.d0
    vg(:,:,:) = 0.d0
    wg(:,:,:) = 0.d0
    !
    if(myid.eq.0) print*, '*** Single phase checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
  elseif(restart_mp_isot) then
    !
    if(myid.eq.0) print*, "=============================================================="
    if(myid.eq.0) print*, " I restart by loading the isothermal multiphase conditition   "
    if(myid.eq.0) print*, "=============================================================="
    !
    is_first_vel = .true.
    is_first_tmp = .true.
    is_first_pos = .true.
    !
    ! load files
    !
    action_load = 'r' 
    call load(action_load,trim(datadir)//'fldu.bin'  ,n,u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldv.bin'  ,n,v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldw.bin'  ,n,w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldp.bin'  ,n,p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvof.bin',n,vof(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldug.bin' ,n,ug(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvg.bin' ,n,vg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldwg.bin' ,n,wg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldtmp.bin',n,tmp(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldsca.bin',n,sca(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(datadir)//'scalar.out',pth,dpthdt_n,time,istep)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at = ', time, 'time step = ', istep, '. ***'
    !
#ifdef USE_VOF
    call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
#endif
#ifdef ENERGY
    call inittmp(initmp,n,dli,phi,vof,tmp)
#else
    tmp(:,:,:) = tmp0
#endif
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    ! 
    pth      = pth_0
    dpthdt_n = 0.d0
    dpthdt   = pth*dpthdt_n
    !
#ifdef USE_VOF
    call update_vof(n,dli,-qb,dzc,dzf,vof,nor,cur,curv,d_thinc)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          rho_gas(i,j,k)   = thermo_rhog(pth,tmp(i,j,k),sca(i,j,k))
          mu_gas(i,j,k)    = thermo_mug(tmp(i,j,k))
          kappa_gas(i,j,k) = thermo_kag(tmp(i,j,k))
          d_lg(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca(i,j,k),1)
          !
          rho(i,j,k)   = rho1*vof(  i,j,k)+(1.d0-vof(i,j,k))*rho_gas(  i,j,k)
          mu(i,j,k)    = mu1*vof(   i,j,k)+(1.d0-vof(i,j,k))*mu_gas(   i,j,k)
          kappa(i,j,k) = kappa1*vof(i,j,k)+(1.d0-vof(i,j,k))*kappa_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho_gas  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu_gas   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa_gas)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_lg     )
    !
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
#else
    rho(0:n(1)+1,0:n(2)+1,0:n(3)+1)       = rho2_0
    mu(0:n(1)+1,0:n(2)+1,0:n(3)+1)        = mu2
    kappa(0:n(1)+1,0:n(2)+1,0:n(3)+1)     = kappa2
    rho_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)   = rho2_0
    mu_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)    = mu2
    kappa_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1) = kappa2
    d_lg(0:n(1)+1,0:n(2)+1,0:n(3)+1)      = d_m12
#endif
    !
#ifdef LOW_MACH
    rho0_min = minval(rho(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,rho0_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
    rho0_min = rho2_0 ! gas phase 
#endif
    !
#ifdef VAP_MASS
    !
    !call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
    call cmpt_normals(n,dli,phi,normlx,normly,normlz)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlx)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normly)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlz)
    !
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmple,-1.d0,0.d0) ! Liq. tmp. extension
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmpge,+1.d0,0.d0) ! Gas. tmp. extension
    !
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca &
#ifdef MULT_COMP
   ,scae_liq,1)                             
#else       
              )
#endif
#ifdef MULT_COMP    
    call initsca(inisca,n,dli,dzci,dzfi,cbcsca,bcsca,rho_gas,d_lg, &
                                        pth,0.d0*u,0.d0*v,0.d0*w,phi,tmpge,tmple,sca2,scae_liq,2)    !new 
#endif
    !
#endif
    !
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca)
    !
    if(myid.eq.0) print*, '*** Isothermal multi phase checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
  elseif(restart_mp_comp) then
    !
    if(myid.eq.0) print*, "============================================="
    if(myid.eq.0) print*, " I restart by loading all the fields "
    if(myid.eq.0) print*, "============================================="
    !
    is_first_vel = .true.
    is_first_tmp = .true.
    is_first_pos = .true.
    !
    ! load the fields
    !
    action_load = 'r' 
    call load(action_load,trim(datadir)//'fldu.bin'  ,n,u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldv.bin'  ,n,v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldw.bin'  ,n,w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldp.bin'  ,n,p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvof.bin',n,vof(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldug.bin' ,n,ug(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldvg.bin' ,n,vg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldwg.bin' ,n,wg(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldtmp.bin',n,tmp(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldsca.bin',n,sca(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldsc2.bin',n,sca2(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(datadir)//'fldscl.bin',n,sca_liq(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(datadir)//'scalar.out',pth,dpthdt_n,time,istep)
    !
    if(myid.eq.0) print*, '*** Checkpoint loaded at = ', time, 'time step = ', istep, '. ***'
    !
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca2)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca_liq)
    !
    if(myid.eq.0) print*, '*** Complete multi phase checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
  endif
#ifdef USE_VOF
  deallocate(pos,xc_p,yc_p,zc_p,rc_p)
#endif
  dpthdt = pth*dpthdt_n
  pold(:,:,:) = p(:,:,:)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pold)
  !
  dudtrkold(:,:,:)   = 0.d0
  dvdtrkold(:,:,:)   = 0.d0
  dwdtrkold(:,:,:)   = 0.d0
  dtmpdtrkold(:,:,:) = 0.d0
  !
#ifdef USE_VOF
  !
  ! update vof quantities and thermophysical properties
  !
  call update_vof(n,dli,-qb,dzc,dzf,vof,nor,cur,curv,d_thinc)
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        rho_gas(i,j,k)   = thermo_rhog(pth,tmp(i,j,k),sca(i,j,k))
        mu_gas(i,j,k)    = thermo_mug(tmp(i,j,k))
        kappa_gas(i,j,k) = thermo_kag(tmp(i,j,k))
        d_lg(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca(i,j,k),1)
#ifdef MULT_COMP
        d_l2g(i,j,k)     = thermo_d_lg(pth,tmp(i,j,k),sca2(i,j,k),2)
        d_l12(i,j,k)     = diff_coef_1_12m
#endif
        !
        rho(i,j,k)   = rho1*vof(  i,j,k)+(1.d0-vof(i,j,k))*rho_gas(  i,j,k)
        mu(i,j,k)    = mu1*vof(   i,j,k)+(1.d0-vof(i,j,k))*mu_gas(   i,j,k)
        kappa(i,j,k) = kappa1*vof(i,j,k)+(1.d0-vof(i,j,k))*kappa_gas(i,j,k)
        !
      enddo
    enddo
  enddo
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho_gas  )
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu_gas   )
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa_gas)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_lg     )
#ifdef MUL_COMP
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l2g     )
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l12     )
#endif
  !
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho  )
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
  !
#else
  rho(0:n(1)+1,0:n(2)+1,0:n(3)+1)       = rho2_0
  mu(0:n(1)+1,0:n(2)+1,0:n(3)+1)        = mu2
  kappa(0:n(1)+1,0:n(2)+1,0:n(3)+1)     = kappa2
  rho_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)   = rho2_0
  mu_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)    = mu2
  kappa_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1) = kappa2
  d_lg(0:n(1)+1,0:n(2)+1,0:n(3)+1)      = d_m12
#endif
  !
#ifdef LOW_MACH
  rho0_min = minval(rho(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,rho0_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
  rho0_min = rho2_0 ! gas phase 
#endif
#ifdef VAR_D_LG
  d_lg_max = maxval(d_lg(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,d_lg_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
  d_lg_max = diff_coef_1_12m!d_m12_2
#endif
#ifdef VAR_KAG
  ka_g_max = maxval(kappa_gas(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,ka_g_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
  ka_g_max = kappa2
#endif
#ifdef VAR_MUG
  mu_g_max = maxval(mu_gas(1:n(1),1:n(2),1:n(3)))
  call mpi_allreduce(MPI_IN_PLACE,mu_g_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
  mu_g_max = mu2
#endif
  !
  !call update_property(n,(/mu1  ,mu2  /),vof,mu  )
  !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu )
  !call update_property(n,(/kappa1,kappa2/),vof,kappa)
  !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
  !call update_property(n,(/cp1   ,cp2   /),vof,cpp   )
  !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,cpp)
  !
  ! compute an initial mass-flux
  !
#ifdef VAP_MASS
  !
  call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
  call cmpt_normals(n,dli,phi,normlx,normly,normlz)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlx)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normly)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlz)
  !
  call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmple,-1.d0,0.d0) ! Liq. tmp. extension
  call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmpge,+1.d0,0.d0) ! Gas. tmp. extension
#ifdef MULT_COMP  
  call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,sca_liq,scae_liq,-1.d0,0.d0)     ! check if zero is fine??????????????
#endif 
  ! 
  ! note: --> tmple is used for Y_int (scae);
  !       --> tmpge is used for rhog and d_lg.
  !
  do k=0,n(3)+1
    do j=0,n(2)+1
      do i=0,n(1)+1
      
        scae(i,j,k) = mass_fraction(pth,tmple(i,j,k) & 
#ifdef MULT_COMP
        ,scae_liq(i,j,k),1)
#else
      ! )
#endif
      enddo
    enddo
  enddo                                                  
  
#ifdef MULT_COMP
  do k=0,n(3)+1
    do j=0,n(2)+1
     do i=0,n(1)+1
       sca2e(i,j,k) = mass_fraction(pth,tmple(i,j,k),scae_liq(i,j,k),2)
     enddo
   enddo
  enddo                                               
#endif
  !
  call cmpt_mflux_m3(n,dli,+0,+0,pth,phi,tmpge,tmple,sca,scae,normlx,normly,normlz,g_sn,+1.d0 &
#ifdef MULT_COMP          
        ,scae_liq,1)  
#else          
          )
#endif          
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,g_sn)
  call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,g_sn,mflux,+1.d0,0.d0)                   
  !
#ifdef MULT_COMP  
  call cmpt_mflux_m3(n,dli,+0,+0,pth,phi,tmpge,tmple,sca2,sca2e,normlx,normly,normlz,g_sn2,+1.d0 &
                                                                                      ,scae_liq,2)  
  !        
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,g_sn2)
  call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,g_sn2,mflux2,+1.d0,0.d0)                     
  !
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        !
        !j_n(i,j,k)  = mflux(i,j,k)*(1.d0-scae(i,j,k))
        !j2_n(i,j,k) = mflux2(i,j,k)*(1.d0-sca2e(i,j,k))
        
        j_n(i,j,k)  = mflux(i,j,k)
        j2_n(i,j,k) = mflux2(i,j,k)
        !
      enddo
    enddo
  enddo
  ! 
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,j_n)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,j2_n)

   do k=1,n(3)
     do j=1,n(2)
       do i=1,n(1)
        !
        mflux_tot(i,j,k)    = -(j_n(i,j,k)+j2_n(i,j,k))/(1.d0-(scae(i,j,k)+sca2e(i,j,k)))
        !
        mflux_a(i,j,k)      = -j_n(i,j,k)+mflux_tot(i,j,k)*scae(i,j,k)
        ! 
        mflux2_a(i,j,k)     = -j2_n(i,j,k)+mflux_tot(i,j,k)*sca2e(i,j,k)
        !
      enddo
    enddo
  enddo
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux_tot)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux_a)
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux2_a)
  !call ab2_scal_liq(n,it_ex,dli,dzci,dzfi,cbcsca,bcsca,dt,rho_liq,d_l12,& 
  !                 !u,v,w,normlx,normly,normlz,phi,mflux,mflux1,sca1_liq)
  !                 u,v,w,normlx,normly,normlz,phi,mflux_tot,mflux_a,sca_liq)
  !call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca_liq)
  !call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,sca_liq,scae_liq,-1.d0,0.d0)
#endif
  !
#else
  ! 
  mflux(1:n(1),1:n(2),1:n(3)) = 10.d0!mflux0
  call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux)
  !
#endif
  !
  ! compute the initial mass of the gas + total interfacial mass-flux
  !
#ifndef MULT_COMP  
  call cmpt_sth(n,dli,vof,rho_gas,mflux,mgas,mfxt) ! not nice name
#else  
  call cmpt_sth(n,dli,vof,rho_gas,mflux_tot,mgas,mfxt) ! not nice name
#endif 
  if(myid.eq.0) print*, 'gas_mass =  ', mgas
  !
  ! visualize field
  !
  !write(fldnum,'(i7.7)') istep
  write(fldnum,'(i12.12)') istep
  include 'out2d.h90'
#ifdef MULT_COMP  
  include 'out3d.h90'
#else
  include 'out3d_single.h90'
#endif  
  i_av   = 0
  istep0 = istep ! we ensure the exact output frequency
  !
  ! initial post-processing
  !
#ifdef USE_VOF
#ifndef VAP_MASS
  call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
#endif
  call cmpt_infq(n,dl,num,rdir_i,phi,sca,tmp,tinf_v,yinf_v)
#ifndef MULT_COMP  
  call int_qtn(n,rho_gas,kappa,tmp,tmpge,tmple,vof,mflux,sca,scae,phi,num,tinf_v,yinf_v,pth,ptho,dpthdt,time,mgas,istep)
#else
  call int_qtn(n,rho_gas,kappa,tmp,tmpge,tmple,vof,mflux_tot,mflux,mflux2,mflux_a,mflux2_a,j_n,j2_n & 
          ,sca,scae,phi,num,tinf_v,yinf_v,pth,ptho,dpthdt,time,mgas,istep,scae_liq,1)
#endif
#endif
  !
  ! compute an initial time-step
  !
  call chkdt(n,dl,u,v,w,rho0_min,d_lg_max,mu_g_max,ka_g_max,dtmax,dtmax_o)
  if(var_dt) then
    dt = min(cfl*dtmax,cfl_o*dtmax_o)
  else
    dt = dt_input 
  endif
  dto = dt
  !
  ! initial postprocessing
  !
  call energy_balance(n,dli,+3,+3,dzc,dzf,time,istep,rho,mu,vof,curv,p, &
                      rho_gas,mu_gas,div_th,divg_th,u,v,w)
  call budget(n,dl,-2,u,v,w,p,rho,rho_gas,vof,time,mgas,istep)
  !
  if(myid.eq.0) print*, 'dtmax = ', dtmax, ' dtmax_o = ', dtmax_o, ' dt = ', dt
  !
  if(myid.eq.0) print*, '==========================================================='
  if(myid.eq.0) print*, '*** Dimensionless physical parameters of the simulation ***'
  if(myid.eq.0) print*, ' re  = ', re_c
  if(myid.eq.0) print*, ' pr  = ', pr_c
  if(myid.eq.0) print*, ' we  = ', we_c
  if(myid.eq.0) print*, ' sc  = ', sc_c
  if(myid.eq.0) print*, ' ste = ', ste_c
  if(myid.eq.0) print*, ' lambda_rho0 = ', lambda_rho0_c
  if(myid.eq.0) print*, ' lambda_mu0  = ', lambda_mu0_c
  if(myid.eq.0) print*, ' lambda_cp0  = ', lambda_cp0_c
  if(myid.eq.0) print*, ' lambda_ka0  = ', lambda_ka0_c
  if(myid.eq.0) print*, ' lambda_mm0  = ', lambda_mm0_c
  if(myid.eq.0) print*, ' lambda_dcp0 = ', lambda_dcp0_c
  if(myid.eq.0) print*, ' phi_p1_c    = ', phi_p1_c
  if(myid.eq.0) print*, ' phi_p2_c    = ', phi_p2_c
  if(myid.eq.0) print*, ' phi_p3_c    = ', phi_p3_c
  if(myid.eq.0) print*, ' phi_p4_c    = ', phi_p4_c
  if(myid.eq.0) print*, '==============================================================='
  if(myid.eq.0) print*, '*** Dimensional physical parameters of the simulation ***'
  if(myid.eq.0) print*, ' uref   = ', uref  , 'lref     = ', lref
  if(myid.eq.0) print*, ' rho1   = ', rho1  , 'rho2     = ', rho2_0
  if(myid.eq.0) print*, ' mu1    = ', mu1   , 'mu2      = ', mu2
  if(myid.eq.0) print*, ' cp1    = ', cp1   , 'cp2      = ', cp2
  if(myid.eq.0) print*, ' kappa1 = ', kappa1, 'kappa2   = ', kappa2
  if(myid.eq.0) print*, ' mm1    = ', m11    , 'm2       = ', m2
  if(myid.eq.0) print*, ' dlg    = ', d_m12 , 'delta_cp = ', delta_cp, 'sigmaca =', sigmaca 
  if(myid.eq.0) print*, ' tmp0    = ', tmpg0 , 'tmpl0      = ', tmpl0 
  if(myid.eq.0) print*, ' lheat  = ', lheat
  if(myid.eq.0) print*, '==============================================================='
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  !
  var(:) = 0.0d0
  v_1(:,:,:) = 0.d0
  v_2(:,:,:) = 0.d0
#ifdef TWOD
      var(3) = sum(vof(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)        !total volume of liquid
#else
      var(3) = sum(vof(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
#endif
      call mpi_allreduce(MPI_IN_PLACE,var(3),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

#ifdef TWOD

      !var(2) = var(2)/(1.d0*(1.d0/4.d0)*pi*lref**2.d0)
      var(2) = ((4.d0*var(3)/pi)**(1.d0/2.d0))/lref
#else
      !var(2) = var(2)/(1.d0*(4.d0/3.d0)*pi*lref**3.d0)
      var(2) = ((6.d0*var(3)/pi)**(1.d0/3.d0))/lref
#endif
      var(1) = time*uref/lref


      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)

      if(sca_liq(i,j,k).eq.0.d0) then

              sca2_liq(i,j,k) = 0.d0
      else
              sca2_liq(i,j,k) = 1.d0-sca_liq(i,j,k)              !only for post processing 
      endif

         enddo
        enddo
      enddo


      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
          var(4)        = var(4) + sca_liq(i,j,k)*vof(i,j,k)   !volume of each species
          var(5)        = var(5) + sca2_liq(i,j,k)*vof(i,j,k)
          enddo
        enddo
      enddo
      
      
      !var(4) = sum(v_1(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)
      !var(5) = sum(v_2(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)



      call mpi_allreduce(MPI_IN_PLACE,var(4),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,var(5),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      nondim = var(4)
      var(6) = var(4)/(var(4)+var(5))!/nondim
      var(7) = var(5)/(var(4)+var(5))!/nondim
      !var(4)  = var(4)*rho_sp_1
      !var(5)  = var(5)*rho_sp_2
      !var(6)  = var(4)/(var(4)+var(5))
      !var(7)  = var(5)/(var(4)+var(5))

      call out0d(trim(datadir)//'volume.out',7,var)



  do while(istep.lt.nstep) 
    !
#ifdef TIMING
    dt12 = MPI_WTIME()
#endif
    !
    istep = istep + 1
    time  = time  + dt
    !
    if(myid.eq.0) print*, '========================================================='
    if(myid.eq.0) print*,  'Timestep #', istep, 'Time = ', time
    if(myid.eq.0) print*, '========================================================='
    !
    ! 0. save the fields from previous time-step
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(rho,pp,p,pold,dt,dtold)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
#ifdef FFT_DIRECT
          pp(i,j,k)   = (1.d0+(dt/dto))*p(i,j,k) - (dt/dto)*pold(i,j,k)
#endif
          pold(i,j,k) = p(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#ifdef FFT_DIRECT
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pp)
#endif
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pold)
    !
#ifdef USE_VOF
    !
    ! 1. vof advection + update of all the quantities
    !
    call advvof(n,dli,dt,-qs,-qb,dzc,dzf,ug,vg,wg,vof,nor,cur,curv,d_thinc)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          rho_gas(i,j,k)   = thermo_rhog(pth,tmp(i,j,k),sca(i,j,k))
          mu_gas(i,j,k)    = thermo_mug(tmp(i,j,k))
          kappa_gas(i,j,k) = thermo_kag(tmp(i,j,k))
          d_lg(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca(i,j,k),1)
#ifdef MULT_COMP          
          d_l2g(i,j,k)      = thermo_d_lg(pth,tmp(i,j,k),sca2(i,j,k),2)
          d_l12(i,j,k)      = diff_coef_1_12m
#endif
          !
          rho(i,j,k)   = rho1*vof(  i,j,k)+(1.d0-vof(i,j,k))*rho_gas(  i,j,k)
          mu(i,j,k)    = mu1*vof(   i,j,k)+(1.d0-vof(i,j,k))*mu_gas(   i,j,k)
          kappa(i,j,k) = kappa1*vof(i,j,k)+(1.d0-vof(i,j,k))*kappa_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho_gas  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu_gas   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa_gas)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_lg     )
#ifdef MULT_COMP
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l2g     )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_l12     )
#endif
    !
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,rho  )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
    !
#else
    rho(0:n(1)+1,0:n(2)+1,0:n(3)+1)       = rho2_0
    mu(0:n(1)+1,0:n(2)+1,0:n(3)+1)        = mu2
    kappa(0:n(1)+1,0:n(2)+1,0:n(3)+1)     = kappa2
    rho_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)   = rho2_0
    mu_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1)    = mu2
    kappa_gas(0:n(1)+1,0:n(2)+1,0:n(3)+1) = kappa2
    d_lg(0:n(1)+1,0:n(2)+1,0:n(3)+1)      = d_m12
#endif
    !
#ifdef LOW_MACH
    rho0_min = minval(rho(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,rho0_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
    rho0_min = rho2_0 ! gas phase 
#endif
    !
    !call update_property(n,(/mu1  ,mu2    /),vof,mu   )
    !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mu   )
    !call update_property(n,(/kappa1,kappa2/),vof,kappa)
    !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
    !call update_property(n,(/cp1   ,cp2   /),vof,cpp  )
    !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,cpp  )
    !
#ifdef VAP_MASS
    !
    ! 2. vapor mass-fraction and mflux --> sca^(n+1), mflux^(n+1)  
    !
    call vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,phi)
    call cmpt_normals(n,dli,phi,normlx,normly,normlz)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlx)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normly)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,normlz)
    !
    call ab2_scal1(n,dli,dzci,dzfi,cbcsca,bcsca,dt,rho_gas,d_lg, &
                                          pth,u,v,w,phi,tmpge,tmple,scae,sca &
#ifdef MULT_COMP
    ,scae_liq,1)
#else
     )
#endif
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca)
    !
    call cmpt_mflux_m3(n,dli,+0,+0,pth,phi,tmpge,tmple,sca,scae,normlx,normly,normlz,g_sn,+1.d0 &
#ifdef MULT_COMP            
    ,scae_liq,1)
#else    
     )
#endif
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,g_sn)
    call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,g_sn,mflux,+1.d0,0.d0)
    !
#ifdef MULT_COMP    
    call ab2_scal1(n,dli,dzci,dzfi,cbcsca,bcsca,dt,rho_gas,d_l2g, &
                                          pth,u,v,w,phi,tmpge,tmple,sca2e,sca2,scae_liq,2)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca2)
    !
    call cmpt_mflux_m3(n,dli,+0,+0,pth,phi,tmpge,tmple,sca2,sca2e,normlx,normly,normlz,g_sn2,+1.d0 &
                                                                      ,scae_liq,2)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,g_sn2)
    call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,g_sn2,mflux2,+1.d0,0.d0)     
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          !j_n(i,j,k)  = mflux(i,j,k)*(1.d0-scae(i,j,k))
          !j2_n(i,j,k) = mflux2(i,j,k)*(1.d0-sca2e(i,j,k))
       
          j_n(i,j,k)  = mflux(i,j,k)
          j2_n(i,j,k) = mflux2(i,j,k)
      
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,j_n)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,j2_n)
    !
    ! compute mflux and the partial ones
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          mflux_tot(i,j,k)    = -(j_n(i,j,k)+j2_n(i,j,k))/(1.d0-(scae(i,j,k)+sca2e(i,j,k)))
          mflux_a(i,j,k)      = -j_n(i,j,k)+mflux_tot(i,j,k)*scae(i,j,k)
          mflux2_a(i,j,k)     = -j2_n(i,j,k)+mflux_tot(i,j,k)*sca2e(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux_tot)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux_a)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux2_a)
#endif
    !
#else
    ! 
    mflux(1:n(1),1:n(2),1:n(3)) = mflux0
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,mflux)
    !
#endif
    !
    ! 4. mass specie in the liquid
    !
#ifdef VAP_MASS
#ifdef MULT_COMP
    ! 4a. liquid species # 1  !contribution 1st species in the liquid phase(robin b.c inside)
    !
    call ab2_scal_liq(istep,n,it_ex,dli,dzci,dzfi,cbcsca,bcsca,dt,rho_liq,d_l12,& 
                       u,v,w,normlx,normly,normlz,phi,mflux_tot,mflux_a,sca_liq)
    call boundp(cbcsca,n,bcsca,dl,dzc,dzf,sca_liq)
    call extended(n,dli,dzc,dzf,+0,phi,normlx,normly,normlz,sca_liq,scae_liq,+1.d0,0.d0)
    !call robin_bc('lin', it_ex,n,dli,dzci,dzfi,-1,mflux_tot,mflux_a,rho_liq,d_l12,phi,normlx,normly,normlz,sca_liq,scae_liq)
    
    !
    !
#endif
#endif

#ifdef ENERGY
    !
    ! 2. temperature advection
    !
#ifndef MULT_COMP    
    call ab2_tmp(is_first_tmp,n,dli,dt,dto,vof,mflux,rho_gas,d_lg,kappa,u,v,w,sca,tmp,dtmpdtrkold,dpthdt)
#else
    call ab2_tmp(is_first_tmp,n,dli,dt,dto,vof,mflux_tot,rho_gas,d_lg,kappa,u,v,w,sca,tmp,dtmpdtrkold,dpthdt)
#endif
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    !
#ifdef VAP_MASS
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmple,-1.d0,0.d0) ! Liq. tmp. extension
    call extended(n,dli,dzc,dzf,+2,phi,normlx,normly,normlz,tmp,tmpge,+1.d0,0.d0) ! Gas. tmp. extension
#else
    tmple(1:n(1),1:n(2),1:n(3)) = tmp(1:n(1),1:n(2),1:n(3))
    call boundp(cbctmp,n,bctmp,dl,dzc,dzf,tmple)
    tmpge(1:n(1),1:n(2),1:n(3)) = tmp(1:n(1),1:n(2),1:n(3))
    call boundp(cbctmp,n,bctmp,dl,dzc,dzf,tmpge)
#endif
    !
#else
    ! 
    tmp(:,:,:) = tmp0
    call boundsb(cbctmp,n,bctmp,dl,dzc,dzf,tmp)
    tmple(:,:,:) = tmpl0
    call boundp(cbctmp,n,bctmp,dl,dzc,dzf,tmple)
    tmpge(:,:,:) = tmp0
    call boundp(cbctmp,n,bctmp,dl,dzc,dzf,tmpge)
    !
#endif
    !
    ! note: --> tmple is used for Y_int (scae);
    !       --> tmpge is used for rhog and d_lg.
    !
#ifdef VAP_MASS
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          scae(i,j,k)  = mass_fraction(pth,tmple(i,j,k) &
#ifdef MULT_COMP          
          ,scae_liq(i,j,k),1)
#else                     
          ) 
#endif
        enddo
      enddo
    enddo
#ifdef MULT_COMP
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          sca2e(i,j,k) = mass_fraction(pth,tmple(i,j,k),scae_liq(i,j,k),2) 
        enddo
      enddo
    enddo
#endif
!    
#else
    sca (:,:,:) = 0.d0
    scae(:,:,:) = 0.d0
#ifdef MULT_COMP
    sca2 (:,:,:) = 0.d0
    sca2e(:,:,:) = 0.d0
#endif    
#endif
    !
#ifdef LOW_MACH
    !
    ! 3. th. pressure advection
    !
    ptho = pth
#ifndef MULT_COMP    
    call cmpt_pth(n,dli,dt,vof,kappa,mflux,rho_gas,d_lg,rho1,tmp,tmpge,sca,scae,mgas,mfxt,pth,dpthdt_n,dpthdt)
#else
    call cmpt_pth(n,dli,dt,vof,kappa,mflux_tot,rho_gas,d_lg,rho1,tmp,tmpge,sca,scae,mgas,mfxt,pth,dpthdt_n,dpthdt)
#endif
    if(myid.eq.0) print*, 'pth  =  ', pth
    if(myid.eq.0) print*, 'mgas =  ', mgas
    !
#else
    !
    pth      = pth_0
    dpthdt_n = 0.d0
    dpthdt   = pth*dpthdt_n
    !
#endif
    !
    ! 4. compute div_th
    !
#ifndef MULT_COMP    
    call cmpt_divth(n,dli,ptho,pth,dpthdt_n,vof,mflux,rho1,kappa,rho_gas,d_lg,tmp,tmpge,sca,scae,div_th,divg_th)
#else
call cmpt_divth(n,dli,ptho,pth,dpthdt_n,vof,mflux_tot,rho1,kappa,rho_gas,d_lg,tmp,tmpge,sca,scae,div_th,divg_th)      !!!check if mflux
#endif
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,div_th )
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,divg_th)
    !
    ! 5. velocity prediction 
    !
    call ab2_mom(istep,is_first_vel,n,dli,dzci,dzfi,dt,dto,u,v,w,tmp,rho_gas,pold,pp,curv,vof,mu,rho,rho0_min, &
                                                      dudtrkold,dvdtrkold,dwdtrkold,up,vp,wp,time)
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp)
    !
#ifdef FFT_DIRECT
    !
    call fillps(n,dli,dzci,dzfi,1.d0/dt,vof,div_th,up,vp,wp,rho,pp,rho0_min,p)
    call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,p)
    call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),p)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    call correc(n,qb,dli,dzci,dt,p,pp,up,vp,wp,rho,rho0_min,u,v,w)
    call bounduvw_b(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
    !
#endif
    !
#ifdef FFT_DIRECT
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(p,pold)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = pold(i,j,k) + p(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
#endif
    !
    ! Interface velocity
    !
#ifdef EVAP
    !
    call fillps(n,dli,dzci,dzfi,1.d0/dt,vof,div_th,0.d0*up,0.d0*vp,0.d0*wp,rho,0.d0*pp,rho0_min,aux)
    call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,p)
    call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),aux)
    call boundp(cbcpre,n,bcpre,dl,dzc,dzf,aux)
    call correc(n,qs,dli,dzci,dt,aux,0.0d0*pp,0.0d0*up,0.0d0*vp,0.0d0*wp,rho,rho0_min,ug,vg,wg)
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
#ifndef MULT_COMP
    call inte_vel(n,nor,mflux,u,v,w,ug,vg,wg)
#else
    call inte_vel(n,nor,mflux_tot,u,v,w,ug,vg,wg)
#endif
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
    !
#else
    !
    ug(1:n(1),1:n(2),1:n(3)) = u(1:n(1),1:n(2),1:n(3))
    vg(1:n(1),1:n(2),1:n(3)) = v(1:n(1),1:n(2),1:n(3))
    wg(1:n(1),1:n(2),1:n(3)) = w(1:n(1),1:n(2),1:n(3))
    call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,ug,vg,wg)
    !
#endif

    dto = dt
    if(mod(istep,icheck).eq.0) then
      !
      ! divergence
      !
      var(:) = 0.0d0
      v_1(:,:,:) = 0.d0
      v_2(:,:,:) = 0.d0
#ifdef TWOD
      var(3) = sum(vof(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)        !total volume of liquid
#else
      var(3) = sum(vof(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
#endif
      call mpi_allreduce(MPI_IN_PLACE,var(3),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

#ifdef TWOD

      !var(2) = var(2)/(1.d0*(1.d0/4.d0)*pi*lref**2.d0)
      var(2) = ((4.d0*var(3)/pi)**(1.d0/2.d0))/lref
#else
      !var(2) = var(2)/(1.d0*(4.d0/3.d0)*pi*lref**3.d0)
      var(2) = ((6.d0*var(3)/pi)**(1.d0/3.d0))/lref
#endif
      var(1) = time*uref/lref


      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)

      if(sca_liq(i,j,k).eq.0.d0) then

              sca2_liq(i,j,k) = 0.d0
      else
              sca2_liq(i,j,k) = 1.d0-sca_liq(i,j,k)              !only for post processing 
      endif

         enddo
        enddo
      enddo


      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)

          var(4) = var(4)+sca_liq(i,j,k)*vof(i,j,k)   !volume of each species
          var(5) = var(5)+sca2_liq(i,j,k)*vof(i,j,k)

          enddo
        enddo
      enddo

      !var(4) = sum(v_1(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)
      !var(5) = sum(v_2(n(1)/2,1:n(2),1:n(3)))*dl(2)*dl(3)



      call mpi_allreduce(MPI_IN_PLACE,var(4),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,var(5),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      var(6) = var(4)/(var(4)+var(5))!/nondim
      var(7) = var(5)/(var(4)+var(5))!/nondim 

      !var(4)  = var(4)*rho_sp_1
      !var(5)  = var(5)*rho_sp_2
      !var(6)  = var(4)/(var(4)+var(5))
      !var(7)  = var(5)/(var(4)+var(5))

      call out0d(trim(datadir)//'volume.out',7,var)

      e_div = 0.d0
      e_div_mean = 0.d0
      do k=1,kmax!-1
        do j=1,jmax!-1
          do i=1,imax!-1
            !
            divu = (u(i,j,k)-u(i-1,j,k))*dxi + &
                   (v(i,j,k)-v(i,j-1,k))*dyi + &
                   (w(i,j,k)-w(i,j,k-1))*dzi
            !
            e_div      = max(abs(div_th(i,j,k)-divu),e_div)
            e_div_mean = e_div_mean + abs(div_th(i,j,k)-divu)*dx*dy*dz/lx/ly/lz
            !
          enddo
        enddo
      enddo                                                                     
      call mpi_allreduce(MPI_IN_PLACE,e_div     ,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,e_div_mean,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,e_div_v)
      if(myid.eq.0) print*, 'Satisfaction of divergence conditions Mean = ', e_div_mean, '| Max = ', e_div
      !
      call chkdiv(n,dli,dzci,u,v,w,divtot,divmax)
      !
      if(isnan(divmax)) then    ! in case of NaN
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
#ifdef VAR_D_LG
      d_lg_max = maxval(d_lg(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,d_lg_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
      d_lg_max = diff_coef_1_12m!d_m12_2
#endif
#ifdef VAR_KAG
      ka_g_max = maxval(kappa_gas(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,ka_g_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
      ka_g_max = kappa2
#endif
#ifdef VAR_MUG
      mu_g_max = maxval(mu_gas(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,mu_g_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#else
      mu_g_max = mu2
#endif
      !
      call chkdt(n,dl,u,v,w,rho0_min,d_lg_max,mu_g_max,ka_g_max,dtmax,dtmax_o)
      if(var_dt) then
        dt = min(cfl*dtmax,cfl_o*dtmax_o)
      else
        dt = dt_input 
      endif
      ! 
      if(myid.eq.0) print*, 'dtmax = ', dtmax, ' dtmax_o = ', dtmax_o, ' dt = ', dt
      !
      !if(dtmax.ne.dtmax) then    ! in case of NaN
      if(isnan(dtmax)) then    ! in case of NaN
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      if(isnan(pth)) then    ! in case of NaN
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
    endif
    !
    if(mod((istep-istep0),iout0d).eq.0) then
      !
      ! general post-processing
      !
#ifdef USE_VOF
      call cmpt_infq(n,dl,num,rdir_i,phi,sca,tmp,tinf_v,yinf_v)
#ifndef MULT_COMP      
      call int_qtn(n,rho_gas,kappa,tmp,tmpge,tmple,vof,mflux,,sca,scae,phi,num,tinf_v,yinf_v,pth,ptho,dpthdt,time,mgas,istep)
#else
      call int_qtn(n,rho_gas,kappa,tmp,tmpge,tmple,vof,mflux_tot,mflux,mflux2,mflux_a,mflux2_a,j_n,j2_n,sca &
                                                  ,scae,phi,num,tinf_v,yinf_v,pth,ptho,dpthdt,time,mgas,istep,scae_liq,1)
#endif
      !
      call cmpt_delta(n,dli,vof,delta)
      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,delta)
#ifndef MULT_COMP      
      call energy_th(n,-2,dl,pth,dpthdt,vof,delta,rho,rho_gas,d_lg,kappa,tmp,sca,mflux,istep,time) 
#else
      call energy_th(n,-2,dl,pth,dpthdt,vof,delta,rho,rho_gas,d_lg,kappa,tmp,sca,mflux_tot,istep,time) 
#endif
      call density_cont(n,-2,dl,pth,dpthdt_n,vof,rho,rho_gas,d_lg,kappa,tmp,sca,istep,time) 
#endif
      !
      ! budget
      !
      call energy_balance(n,dli,+3,+3,dzc,dzf,time,istep,rho,mu,vof,curv,p, &
                          rho_gas,mu_gas,div_th,divg_th,u,v,w)
      call budget(n,dl,-2,u,v,w,p,rho,rho_gas,vof,time,mgas,istep)
      !
    endif
    !
    if(mod((istep-istep0),iout0d_ta).eq.0) then
    !
    endif
    !
    write(fldnum,'(i12.12)') istep
    if(mod((istep-istep0),iout2d).eq.0) then
      include 'out2d.h90'
    endif
    if(mod((istep-istep0),iout3d).eq.0) then
#ifdef MULT_COMP     
            include 'out3d.h90'
#else
            include 'out3d_single.h90'


      p_bottom = 0.d0
      w_bottom = 0.d0
      v_left   = 0.d0
      p_left   = 0.d0      

!   do j=0,jmax+1
!    
!      p_bottom =  (p(0,j,0)+p(0,j,1))/2
!      v_bottom =  (v(0,j,0)+v(0,j,1))/2
!      u_bottom =  (u(0,j,0)+u(0,j,1))/2
!      w_bottom =  w(0,j,0)
!
!      print*,'bottom',0,j, p_bottom,w_bottom,v_bottom,u_bottom  
!
! enddo                                                                     



!  do k=0,kmax+1
!
!      v_left   =  (v(0,0,k)+v(0,1,k))/2
!      p_left   =  (p(0,0,k)+p(0,1,k))/2
!     
!      print*,'left',0,k, p_left,v_left 
!  enddo

#endif
    endif
    if(mod((istep-istep0),isave).eq.0) then
      !
      action_load = 'w' 
      inquire(file='data/fldu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldu.bin   data/fldu_old.bin') ! u
      inquire(file='data/fldv.bin', exist=is_data) 
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldv.bin   data/fldv_old.bin') ! v
      inquire(file='data/fldw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldw.bin   data/fldw_old.bin') ! w
      inquire(file='data/fldp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldp.bin   data/fldp_old.bin') ! p
      inquire(file='data/fldvof.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldvof.bin data/fldvof_old.bin') ! vof
      inquire(file='data/fldug.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldug.bin data/fldug_old.bin') ! ug
      inquire(file='data/fldvg.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldvg.bin data/fldvg_old.bin') ! vg
      inquire(file='data/fldwg.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldwg.bin data/fldwg_old.bin') ! wg
      inquire(file='data/fldtmp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldtmp.bin data/fldtmp_old.bin') ! tmp
      inquire(file='data/fldsca.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldsca.bin data/fldsca_old.bin') ! sca
      inquire(file='data/fldsc2.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldsc2.bin data/fldsc2_old.bin') ! sca2
      inquire(file='data/fldscl.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldscl.bin data/fldscl_old.bin') ! sca_liq
      !
      inquire(file='data/scalar.out', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/scalar.out data/scalar_old.out')
      !
      call load(action_load,trim(datadir)//'fldu.bin'  ,n,u(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldv.bin'  ,n,v(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldw.bin'  ,n,w(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldp.bin'  ,n,p(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldvof.bin',n,vof(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldug.bin' ,n,ug(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldvg.bin' ,n,vg(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldwg.bin' ,n,wg(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldtmp.bin',n,tmp(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldsca.bin',n,sca(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldsc2.bin',n,sca2(1:n(1),1:n(2),1:n(3)))
      call load(action_load,trim(datadir)//'fldscl.bin',n,sca_liq(1:n(1),1:n(2),1:n(3)))
      !
      call load_scalar(action_load,trim(datadir)//'scalar.out',pth,dpthdt_n,time,istep)
      !
      if(myid.eq.0) print*, '*** Checkpoint in fld saved at time = ', time, 'time step = ', istep, '. ***'
      !
    endif
    !
  enddo
  !
  ! 8. clear ffts and/or multigrid
  !
#ifndef MULTI_GRID
  call fftend(arrplanp)
#else
  call destroy_solver_mg
#endif
  !
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
  !
end program cans
