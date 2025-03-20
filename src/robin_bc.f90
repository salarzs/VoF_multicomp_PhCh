module mod_robin_bc
  !
  use mpi
  use decomp_2d
  use mod_common_mpi , only: ierr,myId,coord
  use mod_bound      , only: boundp!,updthalo_1p
  !use mod_rks        , only: rk3_ls
  use mod_funcs      , only: heaviside
  use mod_param     
  use mod_output
  !use mod_extrapl_scal   , only: upwd1,momtnn 
  !
  implicit none
  !
  private
  public robin_bc
  !
  contains
  !
  subroutine robin_bc(iistep,type_pol_reconstr,nstep,n,dli,dzci,dzfi,idir,a,c,rhol,d_liq,phi,normx,normy,normz,s, s_extrpl)
    !
    implicit none
    !
    integer, intent(in) :: iistep
    integer             :: nnstep
    character(len=12) :: fldnum
    integer, intent(in) :: nstep,idir
    integer, intent(in), dimension(3) :: n
    character(len=3), intent(in) :: type_pol_reconstr
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in ), dimension(-2:)         :: dzci,dzfi
    real(8), intent(in), dimension(0:,0:,0:) :: a,c,rhol,d_liq
    real(8), intent(in), dimension(0:,0:,0:) :: normx,normy,normz
    real(8), intent(in), dimension(0:,0:,0:) :: s
    real(8), intent(in), dimension(-2:,-2:,-2:) :: phi
    real(8), intent(out), dimension(0:,0:,0:) :: s_extrpl
    !
    integer :: i,j,k
    integer :: qmin,istep,irk,rkord,iord,rdir
    real(8) :: dt,hi,r
    real(8), dimension(3) :: dl
    real(8), dimension(0:1) :: delta
    real(8), allocatable, dimension(:) :: dzc,dzf
    real(8), allocatable, dimension(:,:,:) :: u,v,w
    real(8), allocatable, dimension(:,:,:) :: interm
    real(8), allocatable, dimension(:,:,:) :: s_extrpl_old,dsedt
    !
    real(8) :: yc,zc
    allocate(dzc(-2:n(3)+3), dzf(-2:n(3)+3))
    allocate(interm(1:n(1),1:n(2),1:n(3)))
    allocate(s_extrpl_old(1:n(1),1:n(2),1:n(3)),dsedt(1:n(1),1:n(2),1:n(3)))
    allocate(u(1:n(1),1:n(2),1:n(3)),v(1:n(1),1:n(2),1:n(3)),w(1:n(1),1:n(2),1:n(3)))
    !
    write(fldnum,'(i12.12)') iistep
    qmin = abs(lbound(normx,1))
    dl = dli**(-1)
    dzc =  dzci**(-1)
    dzf =  dzfi**(-1)
    !
#ifdef TWOD
    dt = 0.49*min(1.d0/dli(2),1.d0/dli(3))
#else
    dt = 0.49*min(1.d0/dli(1),1.d0/dli(2),1.d0/dli(3))
#endif
    rkord = 1 ! order of the rk scheme
    s_extrpl(:,:,:) = 0.d0
    rdir = -1.d0
    do iord = 1,0,-1 
      do k = 1,n(3)
      do j = 1,n(2)
      do i = 1,n(1)
        !
        r = phi(i,j,k)
        u(i,j,k) = normx(i,j,k)*heaviside(r*rdir,0.)*rdir
        v(i,j,k) = normy(i,j,k)*heaviside(r*rdir,0.)*rdir
        w(i,j,k) = normz(i,j,k)*heaviside(r*rdir,0.)*rdir
        interm(i,j,k) = s_extrpl(i,j,k)*heaviside(r*rdir,0.)
        !
      enddo
      enddo
      enddo
     select case(iord)
       case(1)
        call cmpt_ngrad(type_pol_reconstr,n,dli,a,c,rhol,d_liq,phi,normx,normy,normz,s, s_extrpl)  ! compute the normal gradient
        s_extrpl(:,:,:) = s_extrpl(:,:,:)*rdir

       case(0)
         s_extrpl(:,:,:) = s(:,:,:)
     end select
     call boundp(cbcext,n,bcext,dl,dzc,dzf,s_extrpl)

      do istep = 1,40 
        !
        s_extrpl_old(:,:,:) = s_extrpl(1:n(1),1:n(2),1:n(3))
        ! 
        do irk=1,rkord
          call upwd1(n,dli,s_extrpl,u,v,w,dsedt)
          dsedt(:,:,:)=dsedt(:,:,:)+interm(:,:,:)
          call rk3_ls(rkord,irk,n,dt,dsedt,s_extrpl_old,s_extrpl)
          call boundp(cbcext,n,bcext,dl,dzc,dzf,s_extrpl) 
        enddo
        !
      enddo
      !
    enddo
    deallocate(dzc,dzf,interm,s_extrpl_old,dsedt,u,v,w)
    !
    return
  end subroutine robin_bc
  !
  subroutine cmpt_ngrad(type_pol_reconstr,n,dli,a,c,rhol,d_liq,phi,normx,normy,normz,s, ngrad)
    !
    ! computes the linear reconstruction of the normal gradient
    !
    implicit none
    !
    integer, intent(in), dimension(3) :: n
    character(len=3), intent(in) :: type_pol_reconstr
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:,0:,0:) :: a,c,rhol,d_liq
    real(8), intent(in), dimension(0:,0:,0:) :: normx,normy,normz
    real(8), intent(in), dimension(0:,0:,0:) :: s
    real(8), intent(in), dimension(-2:,-2:,-2:) :: phi
    real(8), intent(out), dimension(0:,0:,0:) :: ngrad
    !
    integer :: i,j,k
    real(8) :: b
    !

    select case(type_pol_reconstr)
    case('lin')

      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
             if(phi(i,j,k).ge.0.)then ! compute the normal grad only in the desired domain (liquid)
                if(phi(i,j-1,k).lt.0.d0.or.phi(i,j+1,k).lt.0.d0.or. & 
#ifndef TWOD
                   phi(i-1,j,k).lt.0.d0.or.phi(i+1,j,k).lt.0.d0.or. &
#endif
                   phi(i,j,k-1).lt.0.d0.or.phi(i,j,k+1).lt.0.d0) then
                   !     
                !print*, 'Be careful for the value of diffusion'
                   b = rhol(i,j,k)*4.8157e-7!d_liq(i,j,k)!it's suprising but the "b" coefficient does not match the "c" coeffieint
                                                         !i.e. the c=mflux_1 is computed based on the gradient of vapor mass 
                   ngrad(i,j,k) =  (c(i,j,k)-a(i,j,k)*s(i,j,k))/(-a(i,j,k)*phi(i,j,k)+b+1e-14) ! M. chai et al. [JCP] ! same convention ELE           
                   !
                else ! in the other cells
                   !
                   ngrad(i,j,k) = &
#ifndef TWOD
                                  0.5d0*(s(i+1,j,k)-s(i-1,j,k))*dli(1)*normx(i,j,k)+& 
#endif
                                  0.5d0*(s(i,j+1,k)-s(i,j-1,k))*dli(2)*normy(i,j,k)+&
                                  0.5d0*(s(i,j,k+1)-s(i,j,k-1))*dli(3)*normz(i,j,k) 
                   !
                endif
              else
                      ngrad(i,j,k) = 0.d0
            endif
          enddo
        enddo
      enddo
    case default
            if(myid.eq.0) print*, 'ERROR: invalid name for type_pol_reconstr'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to errors in the case file ***'
      if(myid.eq.0) print*, '    check robin_bc.f90'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    return
  end subroutine cmpt_ngrad
  !
  subroutine rk3_ls(iord,irk,n,dt,dphidt,phiold,phi)
    !
    implicit none
    !
    real(8), parameter, dimension(3) :: alpha1_3 = (/1.d0,3.d0/4.d0,1.d0/3.d0/)
    real(8), parameter, dimension(3) :: alpha2_3 = (/0.d0,1.d0/4.d0,2.d0/3.d0/)
    real(8), parameter, dimension(3) :: beta_3   = (/1.d0,1.d0/4.d0,2.d0/3.d0/)
    real(8), parameter, dimension(3) :: alpha1_2 = (/1.d0,1.d0/2.d0,0.d0/1.d0/)
    real(8), parameter, dimension(3) :: alpha2_2 = (/0.d0,1.d0/2.d0,0.d0/1.d0/)
    real(8), parameter, dimension(3) :: beta_2   = (/1.d0,1.d0/2.d0,0.d0/1.d0/)
    real(8), parameter, dimension(3) :: alpha1_1 = (/1.d0,0.00/1.d0,0.d0/1.d0/)
    real(8), parameter, dimension(3) :: alpha2_1 = (/0.d0,0.d0/1.d0,0.d0/1.d0/)
    real(8), parameter, dimension(3) :: beta_1   = (/1.d0,0.d0/1.d0,0.d0/1.d0/)
    !
    integer, intent(in   )                      :: iord,irk
    integer, intent(in   ), dimension(3)        :: n
    real(8), intent(in   )                      :: dt
    real(8), intent(in   ), dimension( :, :, :) :: dphidt,phiold
    real(8), intent(inout), dimension(0:,0:,0:) :: phi
    !
    real(8), dimension(3) :: alpha1,alpha2,beta
    integer :: i,j,k,p
    !
    select case(iord)
    case(3)
      alpha1(:) = alpha1_3(:)
      alpha2(:) = alpha2_3(:)
      beta(:)   = beta_3(:)
    case(2)
      alpha1(:) = alpha1_2(:)
      alpha2(:) = alpha2_2(:)
      beta(:)   = beta_2(:)
    case(1)
      alpha1(:) = alpha1_1(:)
      alpha2(:) = alpha2_1(:)
      beta(:)   = beta_1(:)
    end select
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          phi(i,j,k) = alpha1(irk)*phiold(i,j,k) + alpha2(irk)*phi(i,j,k)+dt*beta(irk)*dphidt(i,j,k)
        enddo
      enddo
    enddo
    !
    return
  end subroutine rk3_ls
  !
  subroutine upwd1(n,dli,phi,ux,uy,uz,dphidt)
    !
    ! Note: ux,vx,wx are defined already at the cell center
    !       and, therefore, the interpolation is not needed.
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: phi ! phi=te
    real(8), intent(in ), dimension(  :,  :,  :) :: ux,uy,uz
    real(8), intent(out), dimension(  :,  :,  :) :: dphidt
    !
    real(8) :: dphidx,dphidy,dphidz
    integer :: a,i,j,k
    real(8) :: uxc,uyc,uzc
    !


    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#ifndef TWOD
          uxc = ux(i,j,k)
          a = nint(sign(1.d0,uxc))
          dphidx = a*( phi(i,j,k) - phi(i-a,j,k) )*dli(1)
#endif
          !
          uyc = uy(i,j,k)
          a = nint(sign(1.d0,uyc))
          dphidy = a*( phi(i,j,k) - phi(i,j-a,k) )*dli(2)
          !
          uzc = uz(i,j,k)
          a = nint(sign(1.d0,uzc))
          dphidz = a*( phi(i,j,k) - phi(i,j,k-a) )*dli(3)
          !
          dphidt(i,j,k) = -( &
#ifndef TWOD
                  uxc*dphidx + &
#endif
                  uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine upwd1
  !
end module mod_robin_bc
